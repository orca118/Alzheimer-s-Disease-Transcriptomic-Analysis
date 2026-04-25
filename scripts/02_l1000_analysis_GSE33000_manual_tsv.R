library(dplyr)
library(ggplot2)

top_n <- 20

hdac_keywords <- c(
  "vorinostat", "trichostatin", "belinostat",
  "panobinostat", "romidepsin", "HDAC"
)

process_l1000 <- function(file, label) {

  cat("\n============================\n")
  cat("Processing:", label, "\n")
  cat("============================\n")

  df <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)

  # Sort by score
  df <- df %>% arrange(desc(score))

  write.csv(df, paste0(label, "_top_drugs.csv"), row.names = FALSE)

  # ----------------------------
  # Use correct column
  # ----------------------------
  drug_names <- as.character(df$Perturbation)

  # HDAC labeling
  df$HDAC <- grepl(
    paste(hdac_keywords, collapse="|"),
    drug_names,
    ignore.case = TRUE
  )

  # ----------------------------
  # Fisher test
  # ----------------------------
  top_flag <- rep(FALSE, nrow(df))
  top_flag[1:min(top_n, nrow(df))] <- TRUE

  cont_table <- table(Top = top_flag, HDAC = df$HDAC)

  fisher_res <- fisher.test(cont_table)

  cat("\nFisher test result:\n")
  print(fisher_res)

  writeLines(capture.output(fisher_res),
             paste0(label, "_hdac_enrichment_stats.txt"))

  # ----------------------------
  # Plot: Top drugs
  # ----------------------------
  top_df <- df[1:min(15, nrow(df)), ]

  p1 <- ggplot(top_df,
               aes(x = reorder(Perturbation, score), y = score)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Top Candidate Drugs (", label, ")", sep=""),
         x = "",
         y = "L1000 score") +
    theme_minimal()

  ggsave(paste0(label, "_top_drugs.png"), p1, width = 6, height = 4)

  # ----------------------------
  # Plot: HDAC labeling
  # ----------------------------
  top_df$Type <- ifelse(top_df$HDAC, "HDAC inhibitor", "Other")

  p2 <- ggplot(top_df,
               aes(x = reorder(Perturbation, score), y = score, fill = Type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("HDAC inhibitor" = "orange", "Other" = "steelblue")) +
    labs(title = paste("HDAC-related hits (", label, ")", sep=""),
         x = "",
         y = "L1000 score") +
    theme_minimal()

  ggsave(paste0(label, "_hdac_enrichment.png"), p2, width = 6, height = 4)
}

# Run
process_l1000("AD_L1000CDS.tsv", "AD")
process_l1000("HD_L1000CDS.tsv", "HD")