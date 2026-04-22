# Environment Setup Guide

This guide covers all prerequisites for running the **scvi-multiome** pipeline,
from raw FASTQ files through to per-cell-type BigWig files.

---

## 1. Required Conda Environments

The pipeline uses **four existing environments**.  No new environment files are
needed.

| Name | Primary tools | Used by |
|------|--------------|---------|
| `pytorch` | scanpy, scvi-tools, muon, anndata, pysam | Python-based data processing |
| `R` | Seurat, Signac, SeuratDisk | R-based QC steps |
| `base-omics` | deeptools, samtools, bedtools | BigWig generation |
| `snakemake_env` | snakemake, mamba | Workflow orchestration |

### 1.1 Verify `pytorch` environment (add scVI packages if needed)

```bash
conda activate pytorch
python -c "import scanpy, scvi, muon, anndata, pysam; print('OK')"
```

If the above fails, install the missing packages:

```bash
mamba install -c conda-forge -c bioconda \
  scanpy scvi-tools muon anndata pysam
mamba install -c conda-forge \
  scikit-learn scipy pandas matplotlib seaborn
```

### 1.2 Verify `R` environment

```bash
conda activate R
R -e "library(Seurat); library(Signac); library(SeuratDisk)"
```

### 1.3 Verify `base-omics` environment

```bash
conda activate base-omics
bamCoverage --version
samtools --version | head -1
```

### 1.4 Verify `snakemake_env`

```bash
conda activate snakemake_env
snakemake --version
```

---

## 2. Cell Ranger Installation

Cell Ranger is required for the initial FASTQ → count matrix step.

### 2.1 Cell Ranger (RNA)

Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads

```bash
# Verify installation
cellranger --version
# Expected: cellranger 10.0.0 (or later)
```

### 2.2 Cell Ranger ATAC

Download from: https://www.10xgenomics.com/support/software/cell-ranger-atac/downloads

```bash
# Verify installation
cellranger-atac --version
# Expected: cellranger-atac 2.2.0 (or later)
```

Add the unpacked directories to your `PATH` in `~/.bashrc`:

```bash
export PATH=/path/to/cellranger-10.0.0:$PATH
export PATH=/path/to/cellranger-atac-2.2.0:$PATH
```

---

## 3. Reference Genome Files

The `generate_bigwig` workflow requires three files per genome.

### 3.1 Required files checklist

For **mm10** (GRCm38):

| File | Default path in config | Purpose |
|------|------------------------|---------|
| `mm10.fa` | `…/database/mm10/mm10.fa` | Reference FASTA |
| `mm10.chrom.sizes` | `…/database/mm10/mm10.chrom.sizes` | Chromosome lengths |
| `gencode.vM25.annotation.gtf` | `…/database/mm10/gencode.vM25…gtf` | Gene annotation |

### 3.2 How to create missing files

**FASTA** — download from UCSC or GENCODE:

```bash
# UCSC mm10
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
samtools faidx mm10.fa
```

**Chrom sizes** — derived from the FASTA index:

```bash
# From an existing .fai index
cut -f1,2 mm10.fa.fai > mm10.chrom.sizes

# Or directly with samtools
samtools faidx mm10.fa
cut -f1,2 mm10.fa.fai > mm10.chrom.sizes
```

**GTF** — download from GENCODE (mouse vM25, equivalent to GRCm38):

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz
```

### 3.3 Verify paths in config

Run the following check before executing the workflow:

```bash
python - <<'EOF'
import yaml, os
cfg = yaml.safe_load(open("config/config.yaml"))
genome = cfg["reference_genome"]
ref = cfg["genomes"][genome]
all_ok = True
for key, path in ref.items():
    exists = os.path.exists(path)
    status = "✅" if exists else "❌ MISSING"
    print(f"  {key}: {status}")
    print(f"     {path}")
    if not exists:
        all_ok = False
if all_ok:
    print("\nAll reference files found.")
else:
    print("\nSome files are missing — update paths in config/config.yaml")
EOF
```

---

## 4. Update `config/config.yaml`

Open `config/config.yaml` and set:

1. **`reference_genome`** — active genome key (`mm10`, `mm39`, or `hg38`)
2. **`genomes.<name>.fasta`** / **`chrom_sizes`** / **`gtf`** — actual paths on your system
3. **`envs`** — confirm environment names match your conda setup
4. **`wnn_muon_object`** — path to the `.h5mu` produced by `workflows/wnn/`
5. **`rna_anndata`** — path to the RNA `.h5ad`
6. **`atac_fragments_template`** — template pointing to fragment BED files

---

## 5. Verify everything before running

```bash
# Quick pre-flight check
conda activate snakemake_env
snakemake --snakefile workflows/generate_bigwig/snakefile \
          --use-conda --dry-run --cores 1
```

A dry-run without errors means the workflow graph is valid.  The checkpoint
that detects cell types will resolve at runtime once the `.h5mu` file is
available.
