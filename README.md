# scvi-multiome

An end-to-end Snakemake pipeline for processing and integrating single-cell
multiome (RNA + ATAC) data — from raw FASTQ files to per-cell-type BigWig
tracks — using scVI, PoissonVI, and weighted-nearest neighbors (WNN).

## Complete Data Flow

```
Raw FASTQ
  │
  ├─ cellranger count    (RNA)
  ├─ cellranger-atac run (ATAC)
  │
  ▼
workflows/rna/           QC, normalisation, HVG selection
workflows/atac/          QC, peak merging, fragment export
  │
  ▼
workflows/wnn/           scVI + PoissonVI + WNN clustering (Leiden)
  │  → 04_muon_object.h5mu  (cell-type annotations)
  │  → 02_anndata_object_rna.h5ad
  │  → output/{sample}_fragments.bed
  │
  ▼
workflows/generate_bigwig/   per-cell-type BAM → BigWig
  │
  ▼
Results/rna/{celltype}.bw
Results/atac/{celltype}.bw
```

## Key Features

- Automated QC using Gaussian Mixture Models (GMM)
- Separate preprocessing pipelines for RNA and ATAC modalities
- Deep learning feature extraction (scVI / PoissonVI)
- Multimodal integration via WNN
- **Auto-detection** of cell-type annotation column in the `.h5mu` file
- Per-cell-type BigWig generation with `bamCoverage` (RPGC normalisation)
- GPU-accelerated model training (multi-GPU supported)
- Reproducible workflows managed by Snakemake

## Quick Start

See [QUICKSTART.md](QUICKSTART.md) for a 5-step guide.  
See [SETUP.md](SETUP.md) for detailed environment and reference genome setup.

```bash
# 1. Configure
nano config/config.yaml    # set genome paths, env names, input paths
nano config/samples.csv    # add sample IDs and FASTQ locations

# 2. Run upstream workflows
snakemake --snakefile workflows/rna/snakefile  --use-conda --cores 8
snakemake --snakefile workflows/atac/snakefile --use-conda --cores 8
snakemake --snakefile workflows/wnn/snakefile  --use-conda --cores 8

# 3. Generate BigWig files (cell types detected automatically)
snakemake --snakefile workflows/generate_bigwig/snakefile --use-conda --cores 8
```

## Requirements

| Component | Version | Notes |
|-----------|---------|-------|
| Snakemake | ≥ 7 | `snakemake_env` conda env |
| Cell Ranger | ≥ 10 | RNA processing |
| Cell Ranger ATAC | ≥ 2.2 | ATAC processing |
| Python | 3.8+ | `pytorch` conda env |
| scanpy / scvi-tools / muon / anndata / pysam | latest | in `pytorch` env |
| R / Seurat / Signac | latest | `R` conda env |
| deeptools / samtools | latest | `base-omics` conda env |
| GPU | optional | recommended for model training |

## Pipeline Stages

### 1. RNA Processing (`workflows/rna/`)

- Loads raw count matrices → AnnData objects
- GMM-based QC (gene counts, UMIs, MT%, doublets)
- Normalisation, log-transformation, HVG selection
- Merges samples → unified dataset

**Output**: `workflows/rna/objects/…`

### 2. ATAC Processing (`workflows/atac/`)

- Creates ChromatinAssay objects from fragment files
- ATAC-specific QC (TSS enrichment, nucleosome signal, FRiP)
- Merges peaks across samples
- Exports per-sample fragment BED files

**Output**: `workflows/atac/output/{sample_id}_fragments.bed`

### 3. WNN Integration (`workflows/wnn/`)

- Trains scVI (RNA) and PoissonVI (ATAC) models
- Hyperparameter optimisation via Ray Tune
- Computes WNN graph and Leiden clustering
- UMAP visualisation

**Output**: `workflows/wnn/objects/04_muon_object.h5mu`

### 4. BigWig Generation (`workflows/generate_bigwig/`)

- **Auto-detects** cell-type column from the `.h5mu` object  
  (supports `cell_type`, `celltype`, `annotation`, `leiden_celltype`, `leiden`)
- Extracts RNA and ATAC barcodes per cell type
- Converts RNA count matrix → coordinate-sorted BAM
- Converts ATAC fragment records → coordinate-sorted BAM
- Runs `bamCoverage` (deeptools) to produce RPGC-normalised BigWig files

**Output**: `Results/rna/{celltype}.bw`, `Results/atac/{celltype}.bw`

## Configuration

All pipeline parameters live in `config/`:

| File | Purpose |
|------|---------|
| `config/config.yaml` | Reference genomes, conda env names, input/output paths, BigWig params |
| `config/samples.csv` | Single source of truth for sample metadata |
| `config/README_CONFIG.md` | Detailed documentation for every config option |

## Project Structure

```
scvi-multiome/
├── config/
│   ├── config.yaml           # Main configuration
│   ├── samples.csv           # Sample metadata
│   └── README_CONFIG.md      # Config documentation
├── workflows/
│   ├── rna/                  # RNA QC & processing
│   ├── atac/                 # ATAC QC & processing
│   ├── wnn/                  # Multimodal integration
│   ├── generate_bigwig/      # BigWig generation
│   │   ├── snakefile
│   │   └── scripts/
│   │       ├── get_celltypes.py
│   │       ├── extract_cells_by_type.py
│   │       ├── counts_to_bam.py
│   │       └── fragments_to_bam.py
│   ├── celloracle/           # (Future: GRN inference)
│   └── topics/               # (Future: topic modelling)
├── Results/
│   ├── rna/                  # Per-cell-type RNA BigWig files
│   └── atac/                 # Per-cell-type ATAC BigWig files
├── SETUP.md
├── QUICKSTART.md
└── README.md
```

## Methods

### scVI

Scalable probabilistic framework for scRNA-seq using variational autoencoders.  
**Paper**: [Lopez et al., Nature Methods 2018](https://www.nature.com/articles/s41592-018-0229-2)

### PoissonVI

Extension of scVI for scATAC-seq (Poisson distribution for sparse peaks).  
**Paper**: [Martens et al., Nature Methods 2024](https://www.nature.com/articles/s41592-023-02112-6)

### WNN

Multimodal integration via weighted nearest neighbours, implemented in muon.  
**Paper**: [Hao et al., Cell 2021](https://www.sciencedirect.com/science/article/pii/S0092867421005833)

## Citation

If you use this pipeline, please cite the relevant papers for scVI, PoissonVI,
and WNN listed above.

## License

See [LICENSE](LICENSE) for details.
