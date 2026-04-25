# ===============================
# 03_compare_GSE261050_vs_GSE33000_fixed_paths.R
# Primary: GSE261050
# Validation: GSE33000 AD
# Exploratory: GSE33000 HD
# ===============================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# ----------------------------
# Fixed paths from your folders
# ----------------------------
primary_dir <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/data/GSE261050/tables"
valid_dir   <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/GSE33000_l1000_manual_results/tables"
out_dir     <- "/vf/users/leihaiyan/R/alzheimers_drug_repurposing/comparison_results"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
clean_term <- function(x) {
  x <- tolower(as.character(x))
  x <- str_replace_all(x, "\\s*\\(go:\\d+\\)", "")
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_squish(x)
}

first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("Missing required file. Tried:\n", paste(paths, collapse = "\n"))
  ok[1]
}

standardize_drugs <- function(df) {
  if (!("drug" %in% names(df))) {
    if ("Perturbation" %in% names(df)) {
      df$drug <- df$Perturbation
    } else if ("pert_desc" %in% names(df)) {
      df$drug <- df$pert_desc
    } else {
      stop("No usable drug column found. Available columns:\n", paste(names(df), collapse = ", "))
    }
  }
  if (!("score" %in% names(df))) stop("No score column found. Available columns:\n", paste(names(df), collapse = ", "))
  df$score <- suppressWarnings(as.numeric(df$score))
  df
}

get_term_col <- function(df) {
  cand <- c("Term", "Description", "Pathway", "Term name")
  found <- cand[cand %in% names(df)][1]
  if (is.na(found)) stop("Could not find pathway term column. Available columns:\n", paste(names(df), collapse = ", "))
  found
}

top_overlap <- function(df1, df2, n = 10) {
  c1 <- get_term_col(df1)
  c2 <- get_term_col(df2)

  t1 <- df1 |> slice_head(n = n) |> pull(all_of(c1)) |> clean_term()
  t2 <- df2 |> slice_head(n = n) |> pull(all_of(c2)) |> clean_term()

  inter <- intersect(t1, t2)
  uni <- union(t1, t2)

  tibble(
    top_n = n,
    overlap_count = length(inter),
    jaccard = ifelse(length(uni) == 0, NA_real_, length(inter) / length(uni)),
    overlap_terms = paste(inter, collapse = "; ")
  )
}

hdac_flag <- function(x) {
  keywords <- c("vorinostat", "trichostatin", "belinostat", "panobinostat", "romidepsin", "HDAC")
  grepl(paste(keywords, collapse = "|"), as.character(x), ignore.case = TRUE)
}

hdac_summary <- function(df, top_n = 15) {
  df <- df |> arrange(desc(score))
  top <- df |> slice_head(n = top_n)
  top$HDAC <- hdac_flag(top$drug)

  tibble(
    top_n = top_n,
    hdac_count = sum(top$HDAC, na.rm = TRUE),
    hdac_prop = mean(top$HDAC, na.rm = TRUE),
    top_hdac_drugs = paste(top$drug[top$HDAC], collapse = "; ")
  )
}

# ----------------------------
# Files that actually exist in your folders
# ----------------------------
primary_drugs_file <- file.path(primary_dir, "L1000_cleaned_drugs.csv")
primary_go_file    <- file.path(primary_dir, "GO_Biological_Process_2021_results.csv")
primary_kegg_file  <- file.path(primary_dir, "KEGG_2021_Human_results.csv")

valid_ad_drugs_file <- file.path(valid_dir, "AD_top_drugs.csv")
valid_hd_drugs_file <- file.path(valid_dir, "HD_top_drugs.csv")

valid_ad_go_file   <- file.path(valid_dir, "GSE33000_AD_GO_Biological_Process_2021_results.csv")
valid_ad_kegg_file <- file.path(valid_dir, "GSE33000_AD_KEGG_2021_Human_results.csv")

valid_hd_go_file   <- file.path(valid_dir, "GSE33000_HD_GO_Biological_Process_2021_results.csv")
valid_hd_kegg_file <- file.path(valid_dir, "GSE33000_HD_KEGG_2021_Human_results.csv")

# ----------------------------
# Read data
# ----------------------------
primary_drugs <- read_csv(primary_drugs_file, show_col_types = FALSE) |> standardize_drugs()
valid_ad_drugs <- read_csv(valid_ad_drugs_file, show_col_types = FALSE) |> standardize_drugs()
valid_hd_drugs <- read_csv(valid_hd_drugs_file, show_col_types = FALSE) |> standardize_drugs()

primary_go   <- read_csv(primary_go_file, show_col_types = FALSE)
primary_kegg <- read_csv(primary_kegg_file, show_col_types = FALSE)
valid_ad_go  <- read_csv(valid_ad_go_file, show_col_types = FALSE)
valid_ad_kegg <- read_csv(valid_ad_kegg_file, show_col_types = FALSE)

# ----------------------------
# HDAC summaries
# ----------------------------
hdac_compare <- bind_rows(
  hdac_summary(primary_drugs, 15) |> mutate(dataset = "GSE261050", contrast = "AD vs Control"),
  hdac_summary(valid_ad_drugs, 15) |> mutate(dataset = "GSE33000", contrast = "AD vs Control"),
  hdac_summary(valid_hd_drugs, 15) |> mutate(dataset = "GSE33000", contrast = "HD vs Control"),
  hdac_summary(primary_drugs, 20) |> mutate(dataset = "GSE261050", contrast = "AD vs Control"),
  hdac_summary(valid_ad_drugs, 20) |> mutate(dataset = "GSE33000", contrast = "AD vs Control"),
  hdac_summary(valid_hd_drugs, 20) |> mutate(dataset = "GSE33000", contrast = "HD vs Control")
) |>
  mutate(top_n_label = paste0("Top ", top_n)) |>
  select(dataset, contrast, top_n_label, top_n, hdac_count, hdac_prop, top_hdac_drugs)

write_csv(hdac_compare, file.path(out_dir, "hdac_summary.csv"))

# ----------------------------
# Figure 5: AD only (clean version)
# ----------------------------
hdac_ad <- hdac_compare %>%
  filter(contrast == "AD vs Control")

p_hdac <- ggplot(hdac_ad,
                 aes(x = top_n_label, y = hdac_prop, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Transcriptional regulator drugs in Alzheimer’s disease",
    x = "Top-N drugs",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "Fig5_hdac_comparison_AD_only.png"),
  plot = p_hdac,
  width = 6,
  height = 4,
  dpi = 300
)

# ----------------------------
# Pathway overlap: primary AD vs validation AD
# ----------------------------
kegg_overlap_10 <- top_overlap(primary_kegg, valid_ad_kegg, 10)
kegg_overlap_15 <- top_overlap(primary_kegg, valid_ad_kegg, 15)
go_overlap_10   <- top_overlap(primary_go, valid_ad_go, 10)
go_overlap_15   <- top_overlap(primary_go, valid_ad_go, 15)

pathway_overlap <- bind_rows(
  kegg_overlap_10 |> mutate(database = "KEGG", top_n_label = "Top 10"),
  kegg_overlap_15 |> mutate(database = "KEGG", top_n_label = "Top 15"),
  go_overlap_10   |> mutate(database = "GO BP", top_n_label = "Top 10"),
  go_overlap_15   |> mutate(database = "GO BP", top_n_label = "Top 15")
) |>
  select(database, top_n_label, top_n, overlap_count, jaccard, overlap_terms)

write_csv(pathway_overlap, file.path(out_dir, "pathway_overlap.csv"))

p_overlap <- ggplot(pathway_overlap, aes(x = top_n_label, y = overlap_count, fill = database)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Pathway overlap: GSE261050 AD vs GSE33000 AD",
    x = "Top-N pathways",
    y = "Overlap count"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(out_dir, "Fig4_pathway_overlap_comparison.png"),
  plot = p_overlap,
  width = 8,
  height = 4.8,
  dpi = 300
)

# ----------------------------
# Optional: save HD exploratory note if those files exist
# ----------------------------
hd_note <- c(
  "GSE33000 HD exploratory contrast",
  "--------------------------------",
  "HD is treated as a separate exploratory neurodegenerative comparison, not the validation set for AD.",
  paste0("HD top 15 HDAC proportion: ", round(hdac_summary(valid_hd_drugs, 15)$hdac_prop, 3)),
  paste0("HD top 20 HDAC proportion: ", round(hdac_summary(valid_hd_drugs, 20)$hdac_prop, 3))
)
writeLines(hd_note, file.path(out_dir, "GSE33000_HD_note.txt"))

if (file.exists(valid_hd_go_file)) file.copy(valid_hd_go_file, file.path(out_dir, basename(valid_hd_go_file)), overwrite = TRUE)
if (file.exists(valid_hd_kegg_file)) file.copy(valid_hd_kegg_file, file.path(out_dir, basename(valid_hd_kegg_file)), overwrite = TRUE)

# ----------------------------
# Text report
# ----------------------------
report <- c(
  "Comparison Summary",
  "==================",
  "",
  "Primary discovery dataset: GSE261050 AD vs Control",
  "Validation dataset: GSE33000 AD vs Control",
  "Exploratory comparison: GSE33000 HD vs Control",
  "",
  "HDAC summary:",
  paste(capture.output(print(hdac_compare)), collapse = "\n"),
  "",
  "Pathway overlap summary:",
  paste(capture.output(print(pathway_overlap)), collapse = "\n"),
  "",
  "Interpretation note:",
  "- Compare AD to AD across datasets for validation.",
  "- Treat HD as a separate exploratory contrast.",
  "- Do not use HD as the main validation for the AD claim."
)

writeLines(report, file.path(out_dir, "comparison_report.txt"))

cat("Done. Results saved in:\n", out_dir, "\n")