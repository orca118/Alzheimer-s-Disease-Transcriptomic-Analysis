# ===============================
# 05_make_all_figures_publication_quality.R
# Unified script for Figures 1–5
# Shared styling across all figures
# ===============================

suppressPackageStartupMessages({
  pkgs <- c("readr", "dplyr", "ggplot2", "scales", "stringr", "tibble")
  inst <- rownames(installed.packages())
  to_install <- pkgs[!pkgs %in% inst]
  if (length(to_install) > 0) install.packages(to_install)

  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(tibble)
})

# ----------------------------
# Directories
# ----------------------------
primary_dir <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/data/GSE261050/tables"
valid_dir   <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/GSE33000_l1000_manual_results/tables"
compare_dir <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/comparison_results"
out_dir     <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/figures/publication_unified"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Shared theme & palette
# ----------------------------
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 5)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10.5, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 12, 10, 10)
  )

col_up      <- "#D55E00"
col_down    <- "#0072B2"
col_other   <- "#4E79A7"
col_hdac    <- "#E15759"
col_primary <- "#4E79A7"

clean_term <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s*\\(go:\\d+\\)", "")
  x <- str_replace_all(x, "[^A-Za-z0-9]+", " ")
  str_squish(x)
}

is_hdac <- function(x) {
  x <- tolower(as.character(x))
  grepl("statin|vorinostat|entinostat|belinostat|mocetinostat|panobinostat|romidepsin", x, ignore.case = TRUE)
}

pick_file <- function(paths, label) {
  p <- paths[file.exists(paths)][1]
  if (is.na(p) || length(p) == 0) {
    stop("Could not find ", label, ". Checked:\n", paste(paths, collapse = "\n"))
  }
  p
}

# ============================================================
# Figure 1: Volcano plot (primary)
# ============================================================
volcano_file <- pick_file(
  c(
    file.path(primary_dir, "GSE261050_DESeq2_results_mapped.csv"),
    file.path(primary_dir, "results", "GSE261050_DESeq2_results_mapped.csv"),
    file.path(primary_dir, "tables", "GSE261050_DESeq2_results_mapped.csv")
  ),
  "primary DESeq2 mapped results"
)

volc <- read_csv(volcano_file, show_col_types = FALSE)

if (!all(c("logFC", "adj_pvalue") %in% names(volc))) {
  stop("Volcano file must contain columns: logFC, adj_pvalue")
}

volc <- volc |>
  mutate(
    neglog10_padj = -log10(pmax(adj_pvalue, 1e-300)),
    sig = case_when(
      !is.na(adj_pvalue) & adj_pvalue < 0.05 & logFC > 1  ~ "Up",
      !is.na(adj_pvalue) & adj_pvalue < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p1 <- ggplot(volc, aes(x = logFC, y = neglog10_padj)) +
  geom_point(aes(color = sig), alpha = 0.65, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5, color = "grey30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "grey30") +
  scale_color_manual(values = c("Down" = col_down, "NS" = "grey75", "Up" = col_up)) +
  labs(
    title = "Differential expression in Alzheimer’s disease",
    subtitle = "Primary dataset (GSE261050)",
    x = "log2 fold change (AD - Control)",
    y = "-log10(adjusted p-value)",
    color = NULL
  ) +
  theme_pub +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Figure1_volcano_primary.png"), p1, width = 8.5, height = 6.5, dpi = 300)

# ============================================================
# Figure 2: Top candidate drugs (primary)
# ============================================================
primary_drug_file <- pick_file(
  c(
    file.path(primary_dir, "L1000_cleaned_drugs.csv"),
    file.path(primary_dir, "results", "L1000_cleaned_drugs.csv"),
    file.path(primary_dir, "tables", "L1000_cleaned_drugs.csv")
  ),
  "primary cleaned drug table"
)

drug_primary <- read_csv(primary_drug_file, show_col_types = FALSE)

if (!("drug" %in% names(drug_primary))) {
  if ("Perturbation" %in% names(drug_primary)) {
    drug_primary$drug <- drug_primary$Perturbation
  } else if ("Drug" %in% names(drug_primary)) {
    drug_primary$drug <- drug_primary$Drug
  } else {
    stop("Could not find a drug-name column in primary drug table.")
  }
}

if (!("score" %in% names(drug_primary))) {
  if ("Score" %in% names(drug_primary)) {
    drug_primary$score <- drug_primary$Score
  } else {
    stop("Could not find a score column in primary drug table.")
  }
}

drug_primary <- drug_primary |>
  mutate(
    drug = tolower(as.character(drug)),
    score = as.numeric(score),
    hdac = ifelse(is_hdac(drug), "HDAC inhibitor", "Other")
  ) |>
  filter(!is.na(drug), drug != "", !is.na(score)) |>
  arrange(desc(score))

top15_primary <- drug_primary |> slice_head(n = 15)

p2 <- ggplot(top15_primary, aes(x = reorder(drug, score), y = score, fill = hdac)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("HDAC inhibitor" = col_hdac, "Other" = col_primary)) +
  labs(
    title = "Top candidate drugs",
    subtitle = "Primary dataset (GSE261050)",
    x = NULL,
    y = "L1000CDS2 score",
    fill = NULL
  ) +
  theme_pub

ggsave(file.path(out_dir, "Figure2_top_drugs_primary.png"), p2, width = 8.5, height = 6.5, dpi = 300)

# ============================================================
# Figure 3: KEGG pathways (primary)
# ============================================================
kegg_file <- pick_file(
  c(
    file.path(primary_dir, "KEGG_2021_Human_results.csv"),
    file.path(primary_dir, "results", "KEGG_2021_Human_results.csv"),
    file.path(primary_dir, "tables", "KEGG_2021_Human_results.csv")
  ),
  "primary KEGG results"
)

kegg <- read_csv(kegg_file, show_col_types = FALSE)

if (!all(c("Term", "Adjusted.P.value") %in% names(kegg))) {
  stop("KEGG file must contain columns: Term, Adjusted.P.value")
}

kegg <- kegg |>
  mutate(Term = clean_term(Term)) |>
  arrange(Adjusted.P.value) |>
  slice_head(n = 10)

p3 <- ggplot(kegg, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
  geom_col(fill = col_primary, width = 0.75) +
  coord_flip() +
  labs(
    title = "Enriched KEGG pathways in AD (GSE261050)",
    subtitle = "Primary dataset",
    x = NULL,
    y = "-log10(adjusted p-value)"
  ) +
  theme_pub

ggsave(file.path(out_dir, "Figure3_kegg_primary.png"), p3, width = 8.5, height = 6.5, dpi = 300)

# ============================================================
# Figure 4: Pathway overlap (GSE261050 AD vs GSE33000 AD)
# ============================================================
pathway_overlap_file <- pick_file(
  c(
    file.path(compare_dir, "pathway_overlap.csv"),
    file.path(compare_dir, "Fig4_pathway_overlap_comparison.csv")
  ),
  "pathway overlap summary"
)

pathway_overlap <- read_csv(pathway_overlap_file, show_col_types = FALSE)

if (!all(c("database", "top_n_label", "overlap_count") %in% names(pathway_overlap))) {
  stop("Pathway overlap file must contain columns: database, top_n_label, overlap_count")
}

p4 <- ggplot(pathway_overlap, aes(x = top_n_label, y = overlap_count, fill = database)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Pathway overlap: GSE261050 AD vs GSE33000 AD",
    subtitle = "Validation comparison",
    x = "Top-N pathways",
    y = "Overlap count",
    fill = NULL
  ) +
  theme_pub

ggsave(file.path(out_dir, "Figure4_pathway_overlap_comparison.png"), p4, width = 8.5, height = 5.5, dpi = 300)

# ============================================================
# Figure 5: AD-only drug-class comparison across datasets
# ============================================================
hdac_file <- pick_file(
  c(
    file.path(compare_dir, "hdac_summary.csv"),
    file.path(compare_dir, "Fig5_hdac_summary.csv")
  ),
  "HDAC summary"
)

hdac_compare <- read_csv(hdac_file, show_col_types = FALSE)

if (!all(c("dataset", "contrast", "top_n_label", "hdac_prop") %in% names(hdac_compare))) {
  stop("HDAC summary file must contain columns: dataset, contrast, top_n_label, hdac_prop")
}

hdac_ad <- hdac_compare |> filter(contrast == "AD vs Control")

p5 <- ggplot(hdac_ad, aes(x = top_n_label, y = hdac_prop, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Transcriptional regulator drugs in Alzheimer’s disease",
    subtitle = "AD-only comparison across datasets",
    x = "Top-N drugs",
    y = "Proportion",
    fill = NULL
  ) +
  theme_pub

ggsave(file.path(out_dir, "Figure5_hdac_comparison_AD_only.png"), p5, width = 8.5, height = 5.5, dpi = 300)

message("Done. Unified figures saved in: ", out_dir)
