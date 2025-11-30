# scvi-multiome

An end-to-end Snakemake pipeline for processing and integrating single-cell multiome (RNA + ATAC) data using scVI, PoissonVI, and weighted-nearest neighbors (WNN).

## Overview

This pipeline provides a complete workflow for analyzing single-cell multiome datasets, combining gene expression (RNA-seq) and chromatin accessibility (ATAC-seq) data into a unified representation. The workflow leverages state-of-the-art deep learning models and integrates tools from both the R (Seurat/Signac) and Python (Scanpy/scVI) ecosystems.

**Key Features:**
- Automated quality control using Gaussian Mixture Models (GMM)
- Separate preprocessing pipelines for RNA and ATAC modalities
- Deep learning-based feature extraction using scVI and PoissonVI
- Multimodal integration via weighted-nearest neighbors (WNN)
- GPU-accelerated model training
- Reproducible workflows managed by Snakemake

**Prerequisites:**
- Ambient RNA correction (SoupX or CellBender) should be applied before running the pipeline
- Conda/Mamba for environment management
- GPU resources for model training (optional but recommended)

## Pipeline Architecture

The workflow consists of three sequential stages:

### 1. RNA Processing (`workflows/rna/`)

Processes single-cell RNA-seq data through quality control, filtering, and normalization:

- **Preprocessing**: Loads raw count matrices and creates AnnData objects
- **QC Filtering**: Applies GMM-based quality control on metrics including gene counts, UMI counts, mitochondrial/ribosomal percentages, and doublet scores
- **Processing**: Normalizes, transforms, and selects highly variable genes
- **Merging**: Combines filtered samples into a unified dataset

**Output**: Filtered and normalized RNA count matrix ready for deep learning

### 2. ATAC Processing (`workflows/atac/`)

Processes single-cell ATAC-seq data using Signac for chromatin accessibility analysis:

- **Preprocessing**: Creates ChromatinAssay objects from fragment files
- **QC Filtering**: Filters cells based on ATAC-specific metrics (TSS enrichment, nucleosome signal, FRiP scores)
- **Peak Calling**: Merges peaks across samples to create a unified peak set
- **Quantification**: Rebuilds ATAC assay with merged peaks for consistent features
- **Fragment Export**: Generates fragment files for downstream analysis

**Output**: Peak count matrix and fragment files

### 3. WNN Integration (`workflows/wnn/`)

Integrates RNA and ATAC modalities using deep learning and WNN:

- **Feature Filtering**: Selects informative features from each modality
- **Model Training**: 
  - Trains scVI model on RNA data (negative binomial distribution)
  - Trains PoissonVI model on ATAC data (Poisson distribution)
  - Hyperparameter optimization via autotuning
- **Latent Representation**: Extracts low-dimensional embeddings from trained models
- **WNN Integration**: Computes weighted-nearest neighbors across modalities
- **Clustering & Visualization**: Performs Leiden clustering and UMAP projection

**Output**: Integrated multimodal representation with cell clusters and visualizations

## Methods

### scVI (single-cell Variational Inference)

scVI is a scalable probabilistic framework for analyzing scRNA-seq data using variational autoencoders. It provides:
- Batch effect correction across samples
- Dimensionality reduction with biologically interpretable latent spaces
- Efficient processing of large datasets
- Flexible modeling using negative binomial gene expression distributions

**Paper**: [Lopez et al., Nature Methods 2018](https://www.nature.com/articles/s41592-018-0229-2)

### PoissonVI

PoissonVI extends the scVI framework specifically for scATAC-seq data. Unlike scVI, it models sparse, discrete peak accessibility counts using the Poisson distribution, which better captures the characteristics of chromatin accessibility data.

**Paper**: [Martens et al., Nature Methods 2024](https://www.nature.com/articles/s41592-023-02112-6)

### WNN (Weighted Nearest Neighbors)

WNN, implemented via the muon framework, integrates multiple single-cell modalities by computing separate nearest neighbor graphs for each data type and weighting them based on each modality's contribution to cellular identity. This enables:
- Robust multimodal integration
- Improved cell type annotation
- Complementary information capture from gene expression and chromatin accessibility

**Paper**: [Hao et al., Cell 2021](https://www.sciencedirect.com/science/article/pii/S0092867421005833)

## Requirements

- **Snakemake**: Workflow management
- **Conda/Mamba**: Environment management (environments defined in `envs/` directories)
- **Python packages**: scvi-tools, scanpy, muon, pytorch
- **R packages**: Seurat, Signac, Seuratdisk
- **Hardware**: GPU recommended for model training (supports multi-GPU)

## Usage

```bash
# Run the full pipeline (all three stages)
snakemake --use-conda --cores all

# Run individual workflows
snakemake --snakefile workflows/rna/snakefile --use-conda
snakemake --snakefile workflows/atac/snakefile --use-conda
snakemake --snakefile workflows/wnn/snakefile --use-conda
```

## Project Structure

```
scvi-multiome/
├── workflows/
│   ├── rna/              # RNA-seq processing pipeline
│   ├── atac/             # ATAC-seq processing pipeline
│   ├── wnn/              # Multimodal integration pipeline
│   ├── celloracle/       # (Future: gene regulatory network inference)
│   └── topics/           # (Future: topic modeling)
└── README.md
```

## Citation

If you use this pipeline, please cite the relevant papers for scVI, PoissonVI, and WNN listed above.

## License

See LICENSE file for details.
