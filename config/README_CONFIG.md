# Configuration File Documentation

This directory contains the central configuration for the **scvi-multiome** pipeline.  
Edit these files before running any workflow.

---

## `config.yaml` — Main Configuration

### Section 1: `reference_genome`

```yaml
reference_genome: "mm10"
```

The name of the active genome.  Must match one of the keys under `genomes:`.

### Section 2: `genomes`

Each genome entry requires three files:

| Key | Description | How to create if missing |
|-----|-------------|--------------------------|
| `fasta` | Reference FASTA | Download from UCSC / Ensembl / GENCODE |
| `chrom_sizes` | Chromosome sizes file | `samtools faidx genome.fa && cut -f1,2 genome.fa.fai > genome.chrom.sizes` |
| `gtf` | Gene annotation GTF | Download from GENCODE (mouse: vM25+, human: v44+) |

> **Tip**: Use GENCODE annotation rather than UCSC `refGene.gtf` for
> compatibility with `bamCoverage --effectiveGenomeSize`.

### Section 3: `samples_csv`

Path to `config/samples.csv`.  All workflows read from this single file.

### Section 4: `envs`

Conda **environment names** (not `.yml` file paths).

| Key | Default name | Contains |
|-----|-------------|---------|
| `python` | `pytorch` | scanpy, scvi-tools, muon, anndata, pysam |
| `r` | `R` | Seurat, Signac, SeuratDisk |
| `bigwig` | `base-omics` | deeptools, samtools |
| `workflow` | `snakemake_env` | snakemake |

### Section 5: `celltype_column`

The column in `mdata.obs` used to label cell types.

- Set to `null` (default) to **auto-detect** from common names:  
  `cell_type` → `celltype` → `annotation` → `leiden_celltype` → `leiden`
- Set to an explicit name (e.g. `"my_labels"`) to skip auto-detection.

### Section 6: Input paths

| Key | Description |
|-----|-------------|
| `wnn_muon_object` | Path to the final `.h5mu` file from the WNN workflow |
| `rna_anndata` | Path to the RNA `.h5ad` used for model training |
| `atac_fragments_template` | Template for per-sample ATAC fragment BED files |

`{sample_id}` in `atac_fragments_template` is replaced with each ATAC
sample ID found in `samples.csv`.

### Section 7: `bigwig`

| Key | Default | Description |
|-----|---------|-------------|
| `normalization` | `RPGC` | `bamCoverage --normalizeUsing` value |
| `bin_size` | `50` | `--binSize` (bp) |
| `min_mapping_quality` | `20` | `--minMappingQuality` |
| `extra_bamcoverage_params` | `""` | Any additional flags for `bamCoverage` |

---

## `samples.csv` — Sample Table

The **single source of truth** for all sample metadata.

| Column | Description |
|--------|-------------|
| `sample_id` | Unique identifier used in file names |
| `modality` | `rna` or `atac` |
| `organism` | e.g. `mouse`, `human` |
| `stage` | Developmental stage or time-point |
| `fastq_dir` | Directory containing raw FASTQ files |
| `description` | Free-text note (file name hints, etc.) |

### Adding a new sample

1. Add a row to `config/samples.csv`.
2. For RNA samples: place FASTQ files in the directory specified by `fastq_dir`,
   following the naming convention `{sample_id}_f1.fastq.gz` / `{sample_id}_r2.fastq.gz`.
3. For ATAC samples: place FASTQ files matching  
   `{sample_id}_R[123]_001.fastq.gz`.

---

## Checking required files

Run this one-liner to verify all configured paths exist:

```bash
python - <<'EOF'
import yaml, os
cfg = yaml.safe_load(open("config/config.yaml"))
genome = cfg["reference_genome"]
ref = cfg["genomes"][genome]
for key, path in ref.items():
    status = "✅" if os.path.exists(path) else "❌ MISSING"
    print(f"  {key}: {status} — {path}")
EOF
```
