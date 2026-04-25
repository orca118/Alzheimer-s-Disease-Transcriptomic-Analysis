# 🧠 Alzheimer’s Disease Transcriptomic Analysis

This project uses transcriptomic data to identify potential therapeutic targets for Alzheimer’s disease through differential gene expression and drug repurposing analysis.

## 📊 Overview
- Dataset 1: GSE261050 (RNA-seq, primary analysis)
- Dataset 2: GSE33000 (microarray, validation)

### Methods
- Differential expression (DESeq2, limma)
- Pathway enrichment (KEGG)
- Drug repurposing (L1000CDS2)
- Cross-dataset validation

## 📁 Project Structure
- `scripts/` – analysis pipelines  
- `data/` – processed datasets  
- `results/` – output tables and figures  
- `paper/` – manuscript and figures  

## 🚀 How to Run
1. Run primary analysis:
   scripts/01_primary_analysis/01_deseq2_pipeline_GSE261050_v4.R  

2. Run validation:
   scripts/02_validation_analysis/01_limma_pipeline_GSE33000_symbolized_v7.R  

3. Run drug analysis:
   scripts/03_drug_repurposing/02_l1000_analysis_GSE33000_manual_tsv.R  

4. Run comparison:
   scripts/04_comparison/03_compare_GSE261050_vs_GSE33000_v1.R  

5. Generate figures:
   scripts/05_figures/05_make_all_figures_publication_quality.R  

## 📌 Notes
- All datasets are publicly available from GEO  
  https://www.ncbi.nlm.nih.gov/geo/

- L1000CDS2 results were manually downloaded due to API instability

## 📄 License
This project is licensed under the MIT License.