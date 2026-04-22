# Quick Start Guide

Five steps from raw FASTQ to per-cell-type BigWig files.

> **Prerequisite**: All conda environments are installed and reference genome
> files are in place.  See [SETUP.md](SETUP.md) for details.

---

## Step 1 — Configure your project

Edit `config/config.yaml` and `config/samples.csv`.

```bash
# Fill in your sample IDs, FASTQ paths, and reference genome paths
nano config/config.yaml
nano config/samples.csv
```

Key fields in `config.yaml`:

| Field | What to set |
|-------|------------|
| `reference_genome` | `mm10`, `mm39`, or `hg38` |
| `genomes.<name>.chrom_sizes` | Absolute path to `.chrom.sizes` file |
| `genomes.<name>.gtf` | Absolute path to gene annotation GTF |
| `wnn_muon_object` | Path to `04_muon_object.h5mu` from the WNN workflow |
| `rna_anndata` | Path to `02_anndata_object_rna.h5ad` |
| `atac_fragments_template` | Template like `workflows/atac/output/{sample_id}_fragments.bed` |

---

## Step 2 — Run upstream workflows (if not already done)

```bash
conda activate snakemake_env

# RNA processing
snakemake --snakefile workflows/rna/snakefile --use-conda --cores 8

# ATAC processing
snakemake --snakefile workflows/atac/snakefile --use-conda --cores 8

# WNN integration & clustering
snakemake --snakefile workflows/wnn/snakefile --use-conda --cores 8
```

After these complete you should have:
- `workflows/wnn/objects/04_muon_object.h5mu`
- `workflows/wnn/objects/02_anndata_object_rna.h5ad`
- `workflows/atac/output/{sample_id}_fragments.bed`

---

## Step 3 — Inspect cell types in your data

```bash
conda activate pytorch
python -c "
import muon
mdata = muon.read_h5mu('workflows/wnn/objects/04_muon_object.h5mu')
print('obs columns:', list(mdata.obs.columns))
# Look for: cell_type, celltype, annotation, leiden_celltype, leiden
for col in ['cell_type','celltype','annotation','leiden_celltype','leiden']:
    if col in mdata.obs.columns:
        print(f'Found column: {col!r}')
        print(mdata.obs[col].value_counts())
        break
"
```

The `generate_bigwig` workflow **auto-detects** the cell-type column —
no manual configuration needed unless your column has an unusual name.

---

## Step 4 — Generate BigWig files

```bash
conda activate snakemake_env

# Dry-run first (recommended)
snakemake --snakefile workflows/generate_bigwig/snakefile \
          --use-conda --dry-run --cores 1

# Full run
snakemake --snakefile workflows/generate_bigwig/snakefile \
          --use-conda --cores 8
```

The workflow will:
1. Detect cell types from the `.h5mu` file (checkpoint)
2. For each cell type — extract RNA and ATAC barcodes
3. Convert RNA count matrix → sorted BAM
4. Convert ATAC fragment records → sorted BAM
5. Run `bamCoverage` → BigWig (normalised by RPGC)

---

## Step 5 — Inspect outputs

```bash
ls Results/rna/    # *.bw files, one per cell type
ls Results/atac/   # *.bw files, one per cell type
```

Load the `.bw` files in [IGV](https://igv.org) or any genome browser that
supports the BigWig format.

Intermediate BAM files are kept in `workflows/generate_bigwig/intermediate/`
for debugging.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `celltype_column … not found` | Wrong column name | Set `celltype_column` in `config.yaml` |
| `No barcodes found` | Barcode format mismatch | Check if barcodes have a modality prefix (`rna:BC`) |
| `bamCoverage` fails | `base-omics` env missing deeptools | `conda activate base-omics && bamCoverage --version` |
| Empty BigWig | No reads after filtering | Lower `min_mapping_quality` in `config.yaml` |

For more detailed setup instructions see [SETUP.md](SETUP.md).
