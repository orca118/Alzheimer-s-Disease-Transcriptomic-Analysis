#!/usr/bin/env Rscript

# ============================================================
# 01_limma_pipeline_GSE33000_symbolized_v3.R
# Microarray differential expression pipeline for GEO GSE33000
#
# What this script does:
# - downloads/loads GSE33000
# - infers AD, HD, and Control samples from GEO metadata
# - runs limma for AD vs Control and HD vs Control
# - maps probe IDs to gene symbols using GEO feature/platform annotation
# - exports gene-symbol-based ranked lists for downstream L1000/GSEA
#
# Notes:
# - GSE33000 is a microarray dataset, so limma is appropriate.
# - The script avoids manual annotation variables and performs mapping
#   end-to-end inside the pipeline.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ----------------------------
# User settings
# ----------------------------
geo_accession <- "GSE33000"
base_dir <- file.path("data", geo_accession)
raw_dir <- file.path(base_dir, "raw")
tables_dir <- file.path(base_dir, "tables")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

alpha_threshold <- 0.05
lfc_threshold <- 1
n_top_genes <- 250

# ----------------------------
# Helpers
# ----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x), "", x)
  str_to_lower(str_squish(x))
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

extract_sample_text <- function(meta_row) paste(unlist(meta_row), collapse = " | ")

infer_group <- function(txt) {
  txt <- clean_text(txt)
  if (str_detect(txt, regex("huntington|(^|[^a-z])hd([^a-z]|$)|huntington's", ignore_case = TRUE))) return("HD")
  if (str_detect(txt, regex("alzheimer|(^|[^a-z])ad([^a-z]|$)|alzheimer's", ignore_case = TRUE))) return("AD")
  if (str_detect(txt, "control|non-demented|nondemented|normal|healthy|ndhcs")) return("Control")
  NA_character_
}

pick_covariates <- function(meta_df) {
  nms <- names(meta_df)
  age_candidates <- nms[str_detect(nms, regex("age", ignore_case = TRUE))]
  sex_candidates <- nms[str_detect(nms, regex("sex|gender", ignore_case = TRUE))]
  list(
    age_col = if (length(age_candidates) > 0) age_candidates[1] else NA_character_,
    sex_col = if (length(sex_candidates) > 0) sex_candidates[1] else NA_character_
  )
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  x[x %in% c("", "NA", "N/A", "na", "n/a", "null", "NULL", "---", "///")] <- NA_character_
  x
}

split_symbol <- function(x) {
  x <- clean_gene_symbol(x)
  if (is.na(x)) return(NA_character_)
  parts <- unlist(str_split(x, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*|\\s*\\|\\s*"))
  parts <- parts[parts != "" & !is.na(parts)]
  if (length(parts) == 0) return(NA_character_)
  clean_gene_symbol(parts[1])
}

pick_symbol_column <- function(feature_df) {
  nms <- names(feature_df)
  cand <- nms[str_detect(
    nms,
    regex("gene.?symbol|symbol|gene.?name|genename|gene.?title|description|gene.?description|gene assignment", ignore_case = TRUE)
  )]
  if (length(cand) > 0) return(cand[1])
  NA_character_
}

get_platform_table <- function(gpl_id) {
  if (is.na(gpl_id) || !nzchar(gpl_id)) return(NULL)
  gpl_obj <- tryCatch(GEOquery::getGEO(gpl_id), error = function(e) NULL)
  if (is.null(gpl_obj)) return(NULL)
  if (is.list(gpl_obj) && !inherits(gpl_obj, "ExpressionSet")) gpl_obj <- gpl_obj[[1]]
  tbl <- tryCatch(Table(gpl_obj), error = function(e) NULL)
  if (is.null(tbl)) return(NULL)
  as.data.frame(tbl, check.names = FALSE)
}

looks_like_probe_id <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^[0-9]{6,}$|^A_\\d+_P\\d+|^ILMN_", x, perl = TRUE)) > 0.5
}

probe_id_fraction <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_real_)
  mean(grepl("^[0-9]{6,}$|^A_\\d+_P\\d+|^ILMN_", x, perl = TRUE))
}

make_volcano <- function(df, title_text, out_file) {
  p <- ggplot(df, aes(x = logFC, y = neglog10_padj)) +
    geom_point(alpha = 0.30, size = 1) +
    geom_point(data = subset(df, significant), alpha = 0.75, size = 1) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha_threshold), linetype = "dashed") +
    labs(
      title = title_text,
      x = "log2 fold change (case - control)",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal(base_size = 12)
  ggsave(out_file, p, width = 8.5, height = 6.5, dpi = 300)
}

save_top_genes_plot <- function(df, title_text, out_file, top_n = 20, direction = c("Up", "Down")) {
  direction <- match.arg(direction)
  top_df <- if (direction == "Up") {
    df |> filter(!is.na(symbol), logFC > 0) |> slice_head(n = top_n)
  } else {
    df |> filter(!is.na(symbol), logFC < 0) |> slice_head(n = top_n)
  }
  if (nrow(top_df) == 0) return(invisible(NULL))
  p <- ggplot(top_df, aes(x = reorder(symbol, logFC), y = logFC)) +
    geom_col() +
    coord_flip() +
    labs(title = title_text, x = NULL, y = "log2 fold change") +
    theme_minimal(base_size = 12)
  ggsave(out_file, p, width = 8.5, height = 6.5, dpi = 300)
}

run_limma_contrast <- function(expr_mat, meta_df, case_label, control_label = "Control") {
  keep <- meta_df$group %in% c(control_label, case_label)
  meta_sub <- meta_df[keep, , drop = FALSE]
  expr_sub <- expr_mat[, meta_sub$sample_id, drop = FALSE]

  meta_sub$group <- factor(meta_sub$group, levels = c(control_label, case_label))

  covars <- c()
  cov_info <- pick_covariates(meta_sub)

  if (!is.na(cov_info$age_col)) {
    age_vals <- safe_num(meta_sub[[cov_info$age_col]])
    if (sum(!is.na(age_vals)) >= nrow(meta_sub) * 0.75) {
      meta_sub$age <- age_vals
      covars <- c(covars, "age")
    }
  }

  if (!is.na(cov_info$sex_col)) {
    sex_vals <- clean_text(meta_sub[[cov_info$sex_col]])
    sex_vals[sex_vals %in% c("m", "male", "man")] <- "male"
    sex_vals[sex_vals %in% c("f", "female", "woman")] <- "female"
    if (length(unique(na.omit(sex_vals))) >= 2) {
      meta_sub$sex <- factor(sex_vals)
      covars <- c(covars, "sex")
    }
  }

  design_formula <- if (length(covars) > 0) {
    as.formula(paste("~", paste(c(covars, "group"), collapse = " + ")))
  } else {
    ~ group
  }

  design <- model.matrix(design_formula, data = meta_sub)
  fit <- lmFit(expr_sub, design)
  fit2 <- eBayes(fit)

  coef_name <- grep(paste0("group", case_label), colnames(design), value = TRUE)
  if (length(coef_name) == 0) stop("Could not find coefficient for case label: ", case_label)

  topTable(fit2, coef = coef_name[1], number = Inf, sort.by = "P") |>
    rownames_to_column("probe_id")
}

make_gene_level_table <- function(probe_table, feature_map) {
  probe_table |>
    left_join(feature_map, by = "probe_id") |>
    mutate(
      symbol = clean_gene_symbol(symbol),
      symbol = map_chr(symbol, split_symbol),
      neglog10_padj = -log10(pmax(adj.P.Val, 1e-300)),
      significant = !is.na(adj.P.Val) & adj.P.Val < alpha_threshold & abs(logFC) >= lfc_threshold
    ) |>
    filter(!is.na(symbol), symbol != "") |>
    arrange(adj.P.Val, desc(abs(logFC))) |>
    group_by(symbol) |>
    slice_head(n = 1) |>
    ungroup() |>
    arrange(adj.P.Val, desc(abs(logFC)))
}

export_gene_lists <- function(gene_res, case_label) {
  sig <- gene_res |> filter(!is.na(adj.P.Val), adj.P.Val < alpha_threshold)

  up_genes <- sig |>
    filter(logFC > 0) |>
    arrange(adj.P.Val, desc(logFC)) |>
    slice_head(n = n_top_genes) |>
    pull(symbol) |>
    unique() |>
    discard(~ is.na(.x) || .x == "")

  down_genes <- sig |>
    filter(logFC < 0) |>
    arrange(adj.P.Val, logFC) |>
    slice_head(n = n_top_genes) |>
    pull(symbol) |>
    unique() |>
    discard(~ is.na(.x) || .x == "")

  if (looks_like_probe_id(c(up_genes, down_genes))) {
    stop("Exported gene lists still look like probe IDs. The probe-to-symbol mapping failed.")
  }

  writeLines(up_genes, file.path(base_dir, paste0("GSE33000_", case_label, "_up_genes_clean.txt")))
  writeLines(down_genes, file.path(base_dir, paste0("GSE33000_", case_label, "_down_genes_clean.txt")))

  write_csv(
    tibble(
      direction = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
      symbol = c(up_genes, down_genes)
    ),
    file.path(tables_dir, paste0("GSE33000_", case_label, "_symbol_gene_lists.csv"))
  )

  invisible(list(up = up_genes, down = down_genes, sig = sig))
}

resolve_feature_map <- function(gset, expr_mat) {
  feature_df <- fData(gset) |>
    as.data.frame(check.names = FALSE) |>
    rownames_to_column("probe_id")
  feature_df$probe_id <- as.character(feature_df$probe_id)

  message("Feature data columns: ", paste(names(feature_df), collapse = ", "))

  map_entrez_to_symbol <- function(entrez_vec) {
    entrez_vec <- as.character(entrez_vec)
    entrez_vec <- str_extract(entrez_vec, "[0-9]+")
    out <- rep(NA_character_, length(entrez_vec))
    ok <- !is.na(entrez_vec) & nzchar(entrez_vec)
    if (any(ok)) {
      mapped <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(entrez_vec[ok]),
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
      out[ok] <- unname(mapped[entrez_vec[ok]])
    }
    out
  }

  build_map_from_df <- function(df, probe_col, symbol_col = NULL, entrez_col = NULL) {
    out <- data.frame(
      probe_id = as.character(df[[probe_col]]),
      symbol = NA_character_,
      stringsAsFactors = FALSE
    )
    if (!is.null(symbol_col) && symbol_col %in% names(df)) {
      out$symbol <- vapply(df[[symbol_col]], split_symbol, character(1))
    } else if (!is.null(entrez_col) && entrez_col %in% names(df)) {
      out$symbol <- map_entrez_to_symbol(df[[entrez_col]])
    }
    out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
    out$probe_id <- as.character(out$probe_id)
    out$symbol <- as.character(out$symbol)
    unique(out)
  }

  feature_map <- NULL

  symbol_col <- pick_symbol_column(feature_df)
  if (!is.na(symbol_col)) {
    message("Using feature column: ", symbol_col)
    feature_map <- build_map_from_df(feature_df, "probe_id", symbol_col = symbol_col)
  }

  if (is.null(feature_map) || nrow(feature_map) == 0) {
    message("Feature table did not provide usable symbols; trying GPL annotation table ...")
    gpl_id <- annotation(gset)
    message("Platform accession: ", gpl_id)
    annot <- get_platform_table(gpl_id)

    if (!is.null(annot) && nrow(annot) > 0) {
      message("GPL annotation columns: ", paste(names(annot), collapse = ", "))

      if ("ID" %in% names(annot)) {
        annot$probe_id <- as.character(annot[["ID"]])
      } else if (!"probe_id" %in% names(annot)) {
        annot$probe_id <- as.character(rownames(annot))
      } else {
        annot$probe_id <- as.character(annot[["probe_id"]])
      }

      gpl_symbol_col <- pick_symbol_column(annot)
      if (!is.na(gpl_symbol_col)) {
        message("Using GPL column: ", gpl_symbol_col)
        feature_map <- build_map_from_df(annot, "probe_id", symbol_col = gpl_symbol_col)
      } else if ("EntrezGeneID" %in% names(annot)) {
        message("No symbol column found; mapping GPL EntrezGeneID to HGNC symbols ...")
        feature_map <- build_map_from_df(annot, "probe_id", entrez_col = "EntrezGeneID")
      }
    }
  }

  if (is.null(feature_map) || nrow(feature_map) == 0) {
    message("Falling back to probe IDs because no symbol annotation was found.")
    feature_map <- data.frame(
      probe_id = as.character(rownames(expr_mat)),
      symbol = as.character(rownames(expr_mat)),
      stringsAsFactors = FALSE
    )
  }

  feature_map$probe_id <- as.character(feature_map$probe_id)
  feature_map$symbol <- as.character(feature_map$symbol)
  feature_map <- feature_map[!is.na(feature_map$symbol) & nzchar(feature_map$symbol), , drop = FALSE]
  feature_map <- unique(feature_map)

  mapped_rate <- probe_id_fraction(feature_map$symbol)
  if (!is.na(mapped_rate)) {
    message(sprintf("Annotation sanity check: %.1f%% of mapped symbols look like probe IDs.", 100 * mapped_rate))
  }

  if (looks_like_probe_id(feature_map$symbol)) {
    warning("The symbol map still looks like probe IDs. Check the platform annotation table and rerun.")
  }

  feature_map
}
# ----------------------------
# Setup directories
# ----------------------------
ensure_dir(base_dir)
ensure_dir(raw_dir)
ensure_dir(tables_dir)
ensure_dir(results_dir)
ensure_dir(figures_dir)

message("Working directory: ", normalizePath(base_dir, mustWork = FALSE))

# ----------------------------
# 1) Download/load GEO series matrix
# ----------------------------
message("Downloading GEO series matrix for ", geo_accession, " if needed ...")
if (!file.exists(file.path(raw_dir, paste0(geo_accession, "_series_matrix.txt.gz")))) {
  tryCatch(
    GEOquery::getGEOSuppFiles(geo_accession, baseDir = raw_dir, makeDirectory = FALSE),
    error = function(e) message("Supplementary file download skipped: ", conditionMessage(e))
  )
}

gset_list <- GEOquery::getGEO(geo_accession, GSEMatrix = TRUE, getGPL = TRUE)
if (length(gset_list) < 1) stop("Could not load GEO series matrix for ", geo_accession)

gset <- gset_list[[1]]
expr_mat <- Biobase::exprs(gset)
meta <- pData(gset) |> as.data.frame(check.names = FALSE) |> rownames_to_column("sample_id")

meta$sample_text <- apply(meta, 1, extract_sample_text)
meta$group <- vapply(meta$sample_text, infer_group, character(1))
meta <- meta |> filter(!is.na(group))

if (!all(meta$sample_id %in% colnames(expr_mat))) {
  message("Sample IDs do not exactly match expression columns; attempting normalized match.")
  norm <- function(x) str_to_lower(gsub("[^a-z0-9]+", "", x))
  meta$key <- norm(meta$sample_id)
  col_key <- norm(colnames(expr_mat))
  common <- intersect(meta$key, col_key)
  if (length(common) < 6) stop("Could not align metadata with expression columns. Please inspect GEO sample names.")
  meta <- meta[match(common, meta$key), ]
  expr_mat <- expr_mat[, match(common, col_key), drop = FALSE]
  meta$sample_id <- colnames(expr_mat)
} else {
  expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]
}

meta$group <- factor(meta$group, levels = c("Control", "AD", "HD"))
meta$sample_id <- as.character(meta$sample_id)
rownames(meta) <- meta$sample_id

write_csv(meta, file.path(tables_dir, "GSE33000_sample_metadata_inferred.csv"))
writeLines(capture.output(table(meta$group)), file.path(results_dir, "GSE33000_group_counts.txt"))

message("Group counts:")
print(table(meta$group, useNA = "ifany"))
message("Expression matrix: ", nrow(expr_mat), " probes/features x ", ncol(expr_mat), " samples")

# ----------------------------
# 2) Build probe-to-symbol map
# ----------------------------
feature_map <- resolve_feature_map(gset, expr_mat)
write_csv(feature_map, file.path(tables_dir, "GSE33000_probe_symbol_map.csv"))

if (looks_like_probe_id(feature_map$symbol)) {
  message("WARNING: exported symbols still look like probe IDs. Inspect GSE33000_probe_symbol_map.csv before downstream analysis.")
}

# ----------------------------
# 3) Run AD vs Control and HD vs Control
# ----------------------------
comparisons <- c("AD", "HD")
summary_rows <- list()

for (case_label in comparisons) {
  if (sum(meta$group == case_label) < 3) {
    warning("Skipping ", case_label, " because too few case samples were detected.")
    next
  }
  if (sum(meta$group == "Control") < 3) {
    warning("Skipping ", case_label, " because too few control samples were detected.")
    next
  }

  message("Running limma for ", case_label, " vs Control ...")
  probe_res <- run_limma_contrast(expr_mat, meta, case_label = case_label, control_label = "Control")
  gene_res <- make_gene_level_table(probe_res, feature_map)

  out_prefix <- file.path(results_dir, paste0("GSE33000_", case_label, "_vs_Control"))
  table_prefix <- file.path(tables_dir, paste0("GSE33000_", case_label, "_vs_Control"))
  fig_prefix <- file.path(figures_dir, paste0("GSE33000_", case_label, "_vs_Control"))

  write_csv(probe_res, paste0(table_prefix, "_probe_level_results.csv"))
  write_csv(gene_res, paste0(table_prefix, "_gene_level_results.csv"))
  write_csv(gene_res, paste0(out_prefix, "_gene_level_results.csv"))

  volcano_df <- gene_res |>
    mutate(
      neglog10_padj = -log10(pmax(adj.P.Val, 1e-300)),
      significant = !is.na(adj.P.Val) & adj.P.Val < alpha_threshold & abs(logFC) >= lfc_threshold
    )

  make_volcano(
    volcano_df,
    title_text = paste0("Differential expression in GSE33000: ", case_label, " vs Control"),
    out_file = paste0(fig_prefix, "_volcano.png")
  )

  save_top_genes_plot(
    gene_res |> arrange(adj.P.Val, desc(abs(logFC))),
    title_text = paste0("Top upregulated genes: ", case_label, " vs Control"),
    out_file = paste0(fig_prefix, "_top_up_genes.png"),
    top_n = 20,
    direction = "Up"
  )

  save_top_genes_plot(
    gene_res |> arrange(adj.P.Val, abs(logFC)),
    title_text = paste0("Top downregulated genes: ", case_label, " vs Control"),
    out_file = paste0(fig_prefix, "_top_down_genes.png"),
    top_n = 20,
    direction = "Down"
  )

  exported <- export_gene_lists(gene_res, case_label)
  sig <- exported$sig

  summary_rows[[case_label]] <- tibble(
    comparison = paste0(case_label, " vs Control"),
    n_probes_tested = nrow(probe_res),
    n_genes_tested = nrow(gene_res),
    n_sig_genes = nrow(sig),
    n_up = length(exported$up),
    n_down = length(exported$down),
    top_gene = if (nrow(gene_res) > 0) gene_res$symbol[1] else NA_character_
  )
}

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(results_dir, "GSE33000_analysis_summary.csv"))

message("Pipeline completed successfully.")
message("Outputs saved under: ", normalizePath(base_dir, mustWork = FALSE))
