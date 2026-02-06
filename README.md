# Computational Biology Library

A curated collection of tools, pipelines, and frameworks used in:

- Single-cell RNA sequencing (scRNA-seq)
- Bioinformatics analysis
- Machine learning for omics data

This repository serves as:

- A personal knowledge base
- A curated list of essential tools
- A script and pipeline library for computational biology


## Main Categories

- scRNA-seq Core Tools
- Bioinformatics Pipelines
- Machine Learning for Omics


# 1. scRNA-seq Core Tools

## Beginner-friendly tools

These are widely used and ideal starting points.

### Scanpy
- Python framework for scRNA-seq analysis  
- End-to-end pipeline: QC → clustering → DE  
- Link: https://github.com/scverse/scanpy

### Seurat
- Most widely used R package for scRNA-seq  
- Strong integration and visualization tools  
- Link: https://github.com/satijalab/seurat

### Cell Ranger
- Official 10x Genomics processing pipeline  
- FASTQ → gene count matrix  
- Link: https://www.10xgenomics.com/support/software/cell-ranger

### scater
- R package for quality control and visualization  
- Link: https://bioconductor.org/packages/scater

### scran
- Normalization and variance modeling  
- Link: https://bioconductor.org/packages/scran


## Intermediate tools

### Harmony
- Batch correction and dataset integration  
- Works with Seurat and Scanpy  
- Link: https://github.com/immunogenomics/harmony

### BBKNN
- Batch correction using nearest neighbors  
- Designed for Scanpy  
- Link: https://github.com/Teichlab/bbknn

### scVI-tools
- Deep generative models for scRNA-seq  
- Integration, denoising, and annotation  
- Link: https://github.com/scverse/scvi-tools

### Monocle3
- Trajectory and pseudotime analysis  
- Link: https://github.com/cole-trapnell-lab/monocle3

### Slingshot
- Trajectory inference using clustering  
- Link: https://bioconductor.org/packages/slingshot

### PAGA (Scanpy)
- Graph-based trajectory inference  
- Built into Scanpy  
- Link: https://scanpy.readthedocs.io

## Annotation and reference mapping

### SingleR
- Automatic cell type annotation  
- Uses reference datasets  
- Link: https://bioconductor.org/packages/SingleR

### CellTypist
- Fast cell type classification  
- Pretrained immune cell models  
- Link: https://github.com/Teichlab/celltypist

### scmap
- Cell mapping across datasets  
- Link: https://bioconductor.org/packages/scmap


## Doublet detection and QC

### Scrublet
- Doublet detection for scRNA-seq  
- Python-based  
- Link: https://github.com/AllonKleinLab/scrublet

### DoubletFinder
- Doublet detection in Seurat  
- Link: https://github.com/chris-mcginnis-ucsf/DoubletFinder

### SoupX
- Ambient RNA correction  
- Link: https://github.com/constantAmateur/SoupX


## Visualization tools

### cellxgene
- Interactive scRNA-seq visualization  
- Web-based viewer  
- Link: https://github.com/chanzuckerberg/cellxgene

### Loupe Browser
- 10x Genomics visualization tool  
- Link: https://www.10xgenomics.com/products/loupe-browser


# 2. Bioinformatics Pipelines

## General workflow tools

### Snakemake
- Workflow manager for reproducible pipelines  
- Link: https://github.com/snakemake/snakemake

### Nextflow
- Scalable workflow engine  
- Cloud and HPC compatible  
- Link: https://github.com/nextflow-io/nextflow


## Quality control

### FastQC
- Read quality assessment  
- Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

### MultiQC
- Aggregates QC results into one report  
- Link: https://github.com/MultiQC/MultiQC


## Alignment tools

### BWA
- Short read alignment  
- Link: https://github.com/lh3/bwa

### Bowtie2
- Fast alignment for short reads  
- Link: https://github.com/BenLangmead/bowtie2

### STAR
- RNA-seq aligner  
- Link: https://github.com/alexdobin/STAR

### HISAT2
- Fast RNA-seq aligner  
- Link: https://github.com/DaehwanKimLab/hisat2


## Quantification tools

### Salmon
- Alignment-free transcript quantification  
- Link: https://github.com/COMBINE-lab/salmon

### Kallisto
- Fast pseudoalignment-based quantification  
- Link: https://github.com/pachterlab/kallisto

### featureCounts
- Gene-level read counting  
- Link: http://bioinf.wehi.edu.au/featureCounts/


## Variant calling

### GATK
- Industry-standard variant calling toolkit  
- Link: https://github.com/broadinstitute/gatk

### FreeBayes
- Bayesian genetic variant detector  
- Link: https://github.com/freebayes/freebayes


## Metagenomics

### Kraken2
- Taxonomic classification of reads  
- Link: https://github.com/DerrickWood/kraken2

### MetaPhlAn
- Microbial community profiling  
- Link: https://github.com/biobakery/MetaPhlAn

### HUMAnN
- Functional profiling of microbiomes  
- Link: https://github.com/biobakery/humann


# 3. Machine Learning for Omics

## Foundation and deep learning models

### Geneformer
- Transformer model for gene expression  
- Link: https://github.com/broadinstitute/geneformer

### scGPT
- Foundation model for single-cell data  
- Link: https://github.com/bowang-lab/scGPT

### scBERT
- Transformer-based cell representation model  
- Link: https://github.com/TencentAILabHealthcare/scBERT


## Latent variable and probabilistic models

### scVI
- Variational autoencoder for scRNA-seq  
- Integration and denoising  
- Link: https://github.com/scverse/scvi-tools

### scANVI
- Semi-supervised cell annotation  
- Link: https://github.com/scverse/scvi-tools

### totalVI
- Joint RNA and protein modeling  
- Link: https://github.com/scverse/scvi-tools


## Classical ML for omics

### scikit-learn
- General machine learning library  
- Used for clustering, classification, regression  
- Link: https://github.com/scikit-learn/scikit-learn

### XGBoost
- Gradient boosting framework  
- Widely used in genomics prediction tasks  
- Link: https://github.com/dmlc/xgboost


## Deep learning frameworks for biology

### PyTorch
- Deep learning framework used in bio-AI models  
- Link: https://github.com/pytorch/pytorch

### TensorFlow
- Deep learning framework  
- Link: https://github.com/tensorflow/tensorflow


## Graph and network models

### PyTorch Geometric
- Graph neural networks  
- Useful for cell-cell interaction modeling  
- Link: https://github.com/pyg-team/pytorch_geometric

### DGL (Deep Graph Library)
- Graph deep learning framework  
- Link: https://github.com/dmlc/dgl


# Scripts

Reusable scripts and notebooks will be stored in:

```
scripts/
```

Categories:

- `scripts/scrna/`
- `scripts/bioinfo/`
- `scripts/utils/`


# Contributing

To add a new tool or script:

1. Open an issue describing the tool.  
2. Use the templates in the `templates/` folder.  
3. Submit a pull request.


# License

MIT License.

