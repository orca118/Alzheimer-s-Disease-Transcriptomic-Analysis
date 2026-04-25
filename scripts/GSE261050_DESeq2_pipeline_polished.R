#!/usr/bin/env Rscript

# ============================================================
# GSE261050 transcriptomic pipeline (polished)
# 1) Fresh GEO download
# 2) Sample matching + DESeq2 differential expression
# 3) Gene symbol mapping
# 4) Export up/down gene lists for L1000CDS2
# 5) Optional L1000CDS2 / enrichment blocks can be enabled later
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(GEOquery)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ----------------------------
# User settings
# ----------------------------
geo_accession <- "GSE261050"
base_dir <- file.path("data", geo_accession)
raw_dir <- file.path(base_dir, "raw")
extracted_dir <- file.path(base_dir, "extracted")
results_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")

standard_counts_name <- file.path(base_dir, "GSE261050_gene_counts_all.csv.gz")
standard_meta_name <- file.path(base_dir, "GSE261050_metadata.csv.gz")

alpha_threshold <- 0.05
top_per_direction_for_mapping <- 250
top_for_l1000 <- 150

# Set to TRUE later if you want to keep the downstream L1000/API step in this script.
run_l1000_step <- FALSE
run_enrichment_step <- FALSE

# ----------------------------
# Helpers
# ----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

standardize_names <- function(df) {
  names(df) <- names(df) |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("[[:space:]]+", "_") |>
    str_replace_all("[^a-z0-9_]+", "_")
  df
}

normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  x <- gsub("\\.", "_", x)
  x <- gsub("[^a-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

pick_first_existing <- function(patterns, files) {
  for (pat in patterns) {
    hit <- files[stringr::str_detect(basename(files), stringr::regex(pat, ignore_case = TRUE))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

read_any_table <- function(path) {
  b <- basename(path)
  if (stringr::str_detect(b, "\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$")) {
    return(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

coerce_count_matrix <- function(count_df) {
  if (!"gene_id" %in% names(count_df)) {
    names(count_df)[1] <- "gene_id"
  }

  count_df <- count_df |>
    dplyr::mutate(dplyr::across(-gene_id, ~ suppressWarnings(as.numeric(.x))))

  mat <- as.data.frame(count_df)
  rownames(mat) <- mat$gene_id
  mat$gene_id <- NULL
  mat <- as.matrix(mat)

  non_integer <- sum(abs(mat - round(mat)) > 1e-8, na.rm = TRUE)
  if (non_integer > 0) {
    warning("Count matrix contains non-integer values; rounding to nearest integer for DESeq2.")
    mat <- round(mat)
  }

  storage.mode(mat) <- "integer"
  mat
}

map_ensembl_to_symbol <- function(ids) {
  ids_clean <- sub("\\.[0-9]+$", "", ids)

  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ids_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )

  map_df <- dplyr::as_tibble(map_df) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)

  map_df
}

save_plot <- function(plot_obj, path, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(path, plot = plot_obj, width = width, height = height, dpi = dpi)
}

dedupe_keep_order <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- x[nzchar(x)]
  x[!duplicated(x)]
}

# ----------------------------
# Setup directories
# ----------------------------
ensure_dir(base_dir)
ensure_dir(raw_dir)
ensure_dir(extracted_dir)
ensure_dir(results_dir)
ensure_dir(fig_dir)

message("Base directory: ", normalizePath(base_dir, mustWork = FALSE))

# ----------------------------
# 1) Download fresh GEO supplementary files
# ----------------------------
message("Downloading supplementary files from GEO: ", geo_accession)
GEOquery::getGEOSuppFiles(geo_accession, baseDir = raw_dir, makeDirectory = FALSE)

all_files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)
archives <- all_files[stringr::str_detect(tolower(all_files), "\\.(tar|tar\\.gz|tgz)$")]
if (length(archives) > 0) {
  for (a in archives) {
    message("Extracting: ", basename(a))
    untar(a, exdir = extracted_dir)
  }
}

search_pool <- unique(c(all_files, list.files(extracted_dir, recursive = TRUE, full.names = TRUE)))
search_pool <- search_pool[file.exists(search_pool)]

counts_candidate <- pick_first_existing(
  patterns = c("gene_counts.*all.*csv", "counts.*csv", "count.*csv", "matrix.*csv", "counts.*tsv", "count.*tsv"),
  files = search_pool
)
meta_candidate <- pick_first_existing(
  patterns = c("metadata.*csv", "meta.*csv", "sample.*csv", "metadata.*tsv", "meta.*tsv", "sample.*tsv"),
  files = search_pool
)

if (is.na(counts_candidate) || is.na(meta_candidate)) {
  message("Available downloaded files:")
  message(paste(" -", basename(search_pool), collapse = "\n"))
  stop(
    "Could not confidently identify counts and metadata files. Please inspect the downloaded GEO files in: ", raw_dir
  )
}

file.copy(counts_candidate, standard_counts_name, overwrite = TRUE)
file.copy(meta_candidate, standard_meta_name, overwrite = TRUE)
message("Standardized counts file: ", standard_counts_name)
message("Standardized metadata file: ", standard_meta_name)

# ----------------------------
# 2) Load and standardize data
# ----------------------------
counts_raw <- read_any_table(standard_counts_name) |> standardize_names()
meta_raw <- read_any_table(standard_meta_name) |> standardize_names()

if (!"gene_id" %in% names(counts_raw)) {
  names(counts_raw)[1] <- "gene_id"
}

meta <- meta_raw
if ("sample_name" %in% names(meta)) {
  meta$sample_key <- meta$sample_name
} else if ("sample_id" %in% names(meta)) {
  meta$sample_key <- meta$sample_id
} else {
  stop("Metadata must contain sample_name or sample_id.")
}

if ("group" %in% names(meta)) {
  meta <- meta |>
    dplyr::mutate(
      group = tolower(trimws(as.character(group))),
      condition = dplyr::case_when(
        group == "ad" ~ "AD",
        group == "control" ~ "Control",
        TRUE ~ NA_character_
      )
    )
} else if ("condition" %in% names(meta)) {
  meta <- meta |>
    dplyr::mutate(
      condition = dplyr::case_when(
        tolower(trimws(as.character(condition))) == "ad" ~ "AD",
        tolower(trimws(as.character(condition))) == "control" ~ "Control",
        TRUE ~ NA_character_
      )
    )
} else {
  stop("Metadata must contain group or condition.")
}

meta <- meta |>
  dplyr::filter(condition %in% c("AD", "Control")) |>
  dplyr::mutate(sample_key = normalize_sample_id(sample_key)) |>
  dplyr::distinct(sample_key, .keep_all = TRUE)

count_col_names <- names(counts_raw)[-1]
count_key_map <- setNames(count_col_names, normalize_sample_id(count_col_names))
meta_keys <- meta$sample_key
common_keys <- names(count_key_map)[names(count_key_map) %in% meta_keys]

if (length(common_keys) < 4) {
  message("First count keys: ", paste(head(names(count_key_map), 10), collapse = ", "))
  message("First metadata keys: ", paste(head(meta_keys, 10), collapse = ", "))
  stop("Too few matching samples between counts and metadata after filtering.")
}

counts <- counts_raw |>
  dplyr::select(gene_id, dplyr::all_of(unname(count_key_map[common_keys])))

names(counts)[-1] <- common_keys

meta <- meta |>
  dplyr::filter(sample_key %in% common_keys) |>
  dplyr::mutate(sample_key = factor(sample_key, levels = common_keys)) |>
  dplyr::arrange(sample_key) |>
  dplyr::mutate(sample_key = as.character(sample_key))

rownames(meta) <- meta$sample_key
count_mat <- coerce_count_matrix(counts)
count_mat <- count_mat[, meta$sample_key, drop = FALSE]
count_mat <- count_mat[rowSums(count_mat, na.rm = TRUE) > 0, , drop = FALSE]

message("Count matrix dimensions: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")
message("AD samples: ", sum(meta$condition == "AD"), "; Control samples: ", sum(meta$condition == "Control"))

# ----------------------------
# 3) Differential expression with DESeq2
# ----------------------------
meta$condition <- factor(meta$condition, levels = c("Control", "AD"))

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "Control"), alpha = alpha_threshold)

# Create a stable, explicit result table with consistent names
res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene")

if ("log2FoldChange" %in% names(res_df)) {
  res_df$logFC <- res_df$log2FoldChange
} else if ("logFC" %in% names(res_df)) {
  res_df$logFC <- res_df$logFC
} else {
  stop("No fold-change column found in DESeq2 results.")
}

if ("padj" %in% names(res_df)) {
  res_df$adj_pvalue <- res_df$padj
} else if ("adj_pvalue" %in% names(res_df)) {
  res_df$adj_pvalue <- res_df$adj_pvalue
} else {
  stop("No adjusted p-value column found in DESeq2 results.")
}

res_df <- res_df |>
  dplyr::mutate(
    gene = sub("\\.[0-9]+$", "", gene)
  ) |>
  dplyr::select(gene, logFC, pvalue, adj_pvalue, dplyr::everything()) |>
  dplyr::arrange(adj_pvalue, dplyr::desc(abs(logFC)))

readr::write_csv(res_df, file.path(results_dir, "GSE261050_deseq2_results.csv"))
message("Saved DESeq2 results to: ", file.path(results_dir, "GSE261050_deseq2_results.csv"))

volcano_df <- res_df |>
  dplyr::filter(!is.na(logFC), !is.na(adj_pvalue)) |>
  dplyr::mutate(
    neglog10 = -log10(pmax(adj_pvalue, 1e-300)),
    sig = adj_pvalue < alpha_threshold & abs(logFC) >= 1
  )

volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = neglog10)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(data = subset(volcano_df, sig), alpha = 0.7, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha_threshold), linetype = "dashed") +
  labs(
    title = "Differential expression in Alzheimer’s disease",
    x = "log2 fold change (AD - Control)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 12)

save_plot(volcano_plot, file.path(fig_dir, "GSE261050_volcano_plot.png"), width = 8, height = 6)

# ----------------------------
# 4) Build disease signature and map gene IDs to symbols
# ----------------------------
sig <- res_df |>
  dplyr::filter(!is.na(logFC), !is.na(adj_pvalue), adj_pvalue < alpha_threshold)

up_genes_raw <- sig |>
  dplyr::filter(logFC > 0) |>
  dplyr::arrange(adj_pvalue, dplyr::desc(logFC)) |>
  dplyr::slice_head(n = top_per_direction_for_mapping) |>
  dplyr::pull(gene)

dn_genes_raw <- sig |>
  dplyr::filter(logFC < 0) |>
  dplyr::arrange(adj_pvalue, logFC) |>
  dplyr::slice_head(n = top_per_direction_for_mapping) |>
  dplyr::pull(gene)

genes_to_map <- unique(c(up_genes_raw, dn_genes_raw))
map_df <- map_ensembl_to_symbol(genes_to_map)
readr::write_csv(map_df, file.path(results_dir, "GSE261050_gene_symbol_mapping_top.csv"))

res_mapped <- res_df |>
  dplyr::left_join(map_df, by = c("gene" = "ENSEMBL")) |>
  dplyr::mutate(
    symbol = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", toupper(gene), SYMBOL),
    symbol = toupper(trimws(symbol))
  )

sig2 <- res_mapped |>
  dplyr::filter(!is.na(logFC), !is.na(adj_pvalue), adj_pvalue < alpha_threshold)

up <- sig2 |>
  dplyr::filter(logFC > 0) |>
  dplyr::arrange(adj_pvalue, dplyr::desc(logFC)) |>
  dplyr::pull(symbol) |>
  dedupe_keep_order() |>
  head(top_for_l1000)

dn <- sig2 |>
  dplyr::filter(logFC < 0) |>
  dplyr::arrange(adj_pvalue, logFC) |>
  dplyr::pull(symbol) |>
  dedupe_keep_order() |>
  head(top_for_l1000)

if (length(up) < 10 || length(dn) < 10) {
  stop("Too few mapped genes for a reliable L1000CDS2 query.")
}

# ----------------------------
# 5) Export gene lists for L1000CDS2
# ----------------------------
write.table(
  up,
  file = file.path(results_dir, "GSE261050_up_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  dn,
  file = file.path(results_dir, "GSE261050_down_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

signature_df <- tibble::tibble(
  upGenes = paste(up, collapse = ";"),
  downGenes = paste(dn, collapse = ";")
)
readr::write_csv(signature_df, file.path(results_dir, "GSE261050_l1000_signature.csv"))

message("Exported L1000 gene lists:")
message(" - ", file.path(results_dir, "GSE261050_up_genes.txt"))
message(" - ", file.path(results_dir, "GSE261050_down_genes.txt"))

# ----------------------------
# 6) Optional L1000CDS2 step
# ----------------------------
if (run_l1000_step) {
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("To run the L1000 step, install packages httr and jsonlite.")
  }

  query_l1000cds2 <- function(up_genes, dn_genes, tag = "Alzheimer reverse signature") {
    payload <- list(
      data = list(upGenes = up_genes, dnGenes = dn_genes),
      config = list(
        aggravate = FALSE,
        searchMethod = "geneSet",
        share = TRUE,
        combination = FALSE
      ),
      meta = list(
        list(key = "Tag", value = tag),
        list(key = "Cell", value = "bulk brain tissue")
      )
    )

    resp <- httr::POST(
      url = "https://maayanlab.cloud/L1000CDS2/query",
      body = payload,
      encode = "json",
      httr::timeout(120)
    )
    httr::stop_for_status(resp)
    jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  }

  message("Querying L1000CDS2...")
  l1000_raw <- query_l1000cds2(up, dn)
  jsonlite::write_json(l1000_raw, file.path(results_dir, "GSE261050_l1000cds2_response.json"), pretty = TRUE, auto_unbox = TRUE)

  if (!is.null(l1000_raw$topMeta) && length(l1000_raw$topMeta) > 0) {
    ranked <- as_tibble(l1000_raw$topMeta) |>
      dplyr::mutate(rank = dplyr::row_number()) |>
      dplyr::relocate(rank, .before = 1)

    keep_cols <- intersect(
      c("rank", "score", "cell_id", "pert_desc", "pert_id", "pubchem_id", "drugchem_id",
        "pert_time", "pert_time_unit", "pert_dose", "pert_dose_unit", "sig_id", "overlap"),
      names(ranked)
    )
    ranked <- ranked |> dplyr::select(dplyr::all_of(keep_cols))

    readr::write_csv(ranked, file.path(results_dir, "GSE261050_l1000cds2_ranked_drugs.csv"))

    unique_drugs <- ranked |>
      dplyr::arrange(dplyr::desc(score)) |>
      dplyr::group_by(pert_desc) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(dplyr::desc(score))

    readr::write_csv(unique_drugs, file.path(results_dir, "GSE261050_l1000cds2_unique_drugs.csv"))

    top10 <- unique_drugs |>
      dplyr::slice_head(n = 10) |>
      dplyr::mutate(pert_desc = factor(pert_desc, levels = rev(pert_desc)))

    drug_plot <- ggplot(top10, aes(x = pert_desc, y = score)) +
      geom_col() +
      coord_flip() +
      labs(
        title = "Top candidate drugs from L1000CDS2",
        x = NULL,
        y = "L1000CDS2 score"
      ) +
      theme_minimal(base_size = 12)

    save_plot(drug_plot, file.path(fig_dir, "GSE261050_top_drugs_barplot.png"), width = 9, height = 6)
  } else {
    warning("No L1000CDS2 hits returned.")
  }
}

# ----------------------------
# 7) Optional pathway enrichment
# ----------------------------
if (run_enrichment_step) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("To run enrichment, install the enrichR package.")
  }
  dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
  gene_list <- unique(c(up, dn))
  enrich_results <- enrichR::enrichr(gene_list, dbs)

  kegg <- enrich_results[["KEGG_2021_Human"]]
  go <- enrich_results[["GO_Biological_Process_2021"]]

  if (!is.null(kegg)) readr::write_csv(kegg, file.path(results_dir, "GSE261050_kegg_results.csv"))
  if (!is.null(go)) readr::write_csv(go, file.path(results_dir, "GSE261050_go_results.csv"))
}

message("Pipeline complete.")
message("Standardized download files:")
message(" - ", standard_counts_name)
message(" - ", standard_meta_name)
