---
title: Computational Biology Library
---

# Computational Biology Library

A curated, practical library of tools and workflows for:

- **scRNA-seq** (Scanpy/Seurat ecosystem)
- **Bioinformatics pipelines** (QC → alignment → quantification → variants → metagenomics)
- **Machine learning for omics** (foundation models + classical ML)


## Start Here

- **Browse the full library (repo README):**  
  👉 **[Open the main README](../README.md)**

- **Understand scRNA-seq workflows (recommended):**  
  👉 **[scRNA-seq Pipeline Overview](scrna/pipeline-overview.md)**


## Quick Navigation

### 1 scRNA-seq
Core steps: QC → normalization → HVGs → PCA/UMAP → clustering → annotation → integration → DE

- **Pipeline overview:** [docs/scrna/pipeline-overview.md](scrna/pipeline-overview.md)
- (Add more pages later) QC • Integration • DE • Trajectory


### 2 Bioinformatics Pipelines
From raw sequencing reads to biological insight.

Suggested sections you can add later:
- QC reporting (FastQC/MultiQC)
- Alignment + quantification (STAR/HISAT2/Salmon)
- Variant calling (GATK)
- Metagenomics (Kraken2/MetaPhlAn)


### 3 Machine Learning for Omics
Modern ML for gene expression and multi-omics, including transformers and generative models.

Suggested sections you can add later:
- Representation learning
- Reference mapping and annotation models
- Foundation models (Geneformer, scGPT)


## Repository Structure

- **`README.md`** → Full curated list of tools + descriptions  
- **`docs/`** → Practical notes + playbooks (this website)  
- **`scripts/`** → Reusable scripts and helpers  
- **`templates/`** → Templates for adding tools/scripts consistently  


## How to Use This Library

1. Use **README** to discover tools.
2. Use **docs** to follow workflow guidance.
3. Use **scripts** to run quick analyses or reuse code blocks.
4. Add new tools using a consistent format.


## Add a New Tool (quick)
If you find a useful repository/tool:

1. Open an Issue in GitHub (recommended)
2. Add it to the README under the right section
3. Commit and push


## About

Maintained by **Taufia Hussain**  
GitHub: **[TaufiaHussain](https://github.com/TaufiaHussain)**

> Goal: keep a clean, modern library that’s useful for research, teaching, and reproducible analysis.
