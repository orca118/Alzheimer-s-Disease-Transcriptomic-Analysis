# Identifying Potential Drug Candidates for Alzheimer’s Disease Using Transcriptomic Analysis

This repository contains the analysis code used for the study:

**Identifying Potential Drug Candidates for Alzheimer’s Disease Using Transcriptomic Analysis**

**Sophia Wang**  
Richard Montgomery High School, USA

*The National High School Journal of Science*, 2026  
Received: April 30, 2026  
Accepted: July 25, 2026  
Electronic access: September 30, 2026

---

## Overview

This study used publicly available transcriptomic datasets to investigate Alzheimer’s disease (AD)-associated gene-expression changes and identify potential drug-repurposing candidates.

The workflow included:

- Differential gene-expression analysis
- KEGG pathway enrichment
- Computational drug repurposing using L1000CDS2
- Independent cross-dataset validation
- Comparison of pathway- and drug-level signals across transcriptomic platforms

The primary analysis used RNA-seq dataset **GSE261050**, and validation was performed using microarray dataset **GSE33000**.

The study was designed to evaluate the reproducibility of transcriptomic disease signatures and computational drug-repurposing predictions across independent Alzheimer’s disease datasets.

---

## Datasets

### GSE261050 — Primary RNA-seq Dataset

The primary dataset contains postmortem bulk RNA-seq data from the anterior cingulate cortex (BA32) and insula.

For the analysis reported in the paper:

- **101 samples**
- **60 unique donors**
- **76 Alzheimer’s disease samples**
- **25 control samples**

Only samples annotated as Alzheimer’s disease or control were included in the primary differential-expression analysis.

Available metadata included:

- Age
- Sex
- Brain region
- Post-mortem interval (PMI)
- RNA integrity number (RIN)
- Sequencing batch
- Extraction batch

The primary DESeq2 model used disease status as the main variable.

### GSE33000 — Independent Validation Dataset

GSE33000 is a microarray dataset containing samples from:

- Alzheimer’s disease
- Huntington’s disease
- Controls

Only Alzheimer’s disease and control samples were used for the primary validation analysis.

Huntington’s disease samples were analyzed separately as an exploratory comparison and were not included in the main AD validation.

Because GSE261050 is RNA-seq and GSE33000 is microarray-based, the validation focused on cross-platform reproducibility of pathway- and drug-level signals rather than exact replication of individual differential-expression results.

---

## Repository Structure

```text
.
├── scripts/
│   ├── 01_primary_analysis/
│   ├── 02_validation_analysis/
│   ├── 03_drug_repurposing/
│   ├── 04_comparison/
│   └── 05_figures/
├── data/
├── results/
├── paper/
└── README.md
