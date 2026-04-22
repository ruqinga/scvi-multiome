"""
counts_to_bam.py
----------------
Convert a subset of an RNA count matrix (H5AD) into a coordinate-sorted BAM
file suitable for use with ``bamCoverage``.

Strategy
~~~~~~~~
For each cell × gene entry with count > 0 the script emits synthetic
single-end reads.  The number of reads emitted equals the UMI count.
Each read is placed at the transcription-start site (TSS) of the gene
(or at the first base of the gene body when a GTF is not supplied).

Gene coordinates are taken from ``adata.var`` when columns ``chr``,
``start``, ``end`` (and optionally ``strand``) are present, otherwise
from a GTF file configured in ``config/config.yaml``.

Usage (via Snakemake)
~~~~~~~~~~~~~~~~~~~~~
    snakemake.input.anndata      — RNA H5AD (all cells)
    snakemake.input.barcodes     — text file, one barcode per line
    snakemake.output.bam         — output sorted + indexed BAM
    snakemake.params.chrom_sizes — path to .chrom.sizes file
    snakemake.params.gtf         — path to GTF (may be empty string)
    snakemake.params.read_length — synthetic read length (default 100)
"""

import sys
import os
import tempfile

import pysam
import scanpy
import pandas

# SAM flag constants
_FLAG_REVERSE = 0x10  # read mapped on reverse strand


READ_LENGTH_DEFAULT = 100


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_chrom_sizes(path):
    """Return {chrom: length} dict from a .chrom.sizes file."""
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def load_gene_coords_from_gtf(gtf_path, gene_names):
    """Extract gene coordinates from a GTF file for the requested genes.

    Returns a dict: gene_name → (chrom, tss, strand)
    where tss is 0-based.
    """
    gene_set = set(gene_names)
    coords = {}
    print(f"[counts_to_bam] Parsing GTF: {gtf_path}")
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = parts[8]
            name = None
            for field in attr.split(";"):
                field = field.strip()
                if field.startswith("gene_name"):
                    name = field.split('"')[1] if '"' in field else field.split()[-1]
                    break
            if name and name in gene_set and name not in coords:
                chrom = parts[0]
                start = int(parts[3]) - 1  # GTF is 1-based → 0-based
                end = int(parts[4]) - 1
                strand = parts[6]
                tss = start if strand == "+" else end
                coords[name] = (chrom, tss, strand)
    return coords


def build_bam_header(chrom_sizes):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": length} for chrom, length in chrom_sizes.items()],
    }
    return pysam.AlignmentHeader.from_dict(header)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(anndata_path, barcodes_path, output_bam, chrom_sizes_path, gtf_path, read_length):
    print(f"[counts_to_bam] Reading barcodes: {barcodes_path}")
    with open(barcodes_path) as fh:
        barcodes = [line.strip() for line in fh if line.strip()]

    print(f"[counts_to_bam] {len(barcodes)} barcodes selected")

    print(f"[counts_to_bam] Reading AnnData: {anndata_path}")
    adata = scanpy.read_h5ad(anndata_path)

    # Subset to selected barcodes
    common = [bc for bc in barcodes if bc in adata.obs_names]
    if not common:
        # Try stripping modality prefix (e.g. "rna:BC" → "BC")
        strip_map = {bc.split(":")[-1]: bc for bc in adata.obs_names}
        common_stripped = [bc for bc in barcodes if bc.split(":")[-1] in strip_map]
        common = [strip_map[bc.split(":")[-1]] for bc in common_stripped]

    if not common:
        raise ValueError(
            "No barcodes from the barcode list were found in the AnnData object.\n"
            f"First 5 barcodes requested: {barcodes[:5]}\n"
            f"First 5 barcodes in AnnData: {list(adata.obs_names[:5])}"
        )

    adata_sub = adata[common, :]
    print(f"[counts_to_bam] Subset AnnData: {adata_sub.shape}")

    # Gene coordinates
    chrom_sizes = parse_chrom_sizes(chrom_sizes_path)

    if all(c in adata.var.columns for c in ["chr", "start", "end"]):
        print("[counts_to_bam] Using gene coordinates from adata.var")
        var_coords = adata.var[["chr", "start", "end"]].copy()
        if "strand" in adata.var.columns:
            var_coords["strand"] = adata.var["strand"]
        else:
            var_coords["strand"] = "+"
        gene_coords = {
            g: (row["chr"], int(row["start"]), row["strand"])
            for g, row in var_coords.iterrows()
        }
    elif gtf_path:
        gene_coords = load_gene_coords_from_gtf(gtf_path, adata.var_names.tolist())
    else:
        raise ValueError(
            "Gene coordinates not found in adata.var (need 'chr', 'start', 'end' columns) "
            "and no GTF path was provided."
        )

    # Write BAM to a temp file, then sort + index
    tmp_unsorted = output_bam + ".unsorted.bam"
    header = build_bam_header(chrom_sizes)

    print(f"[counts_to_bam] Writing synthetic reads …")
    written = 0
    skipped_genes = 0

    # Densify in chunks to avoid memory issues
    import scipy.sparse as sp
    X = adata_sub.X
    if sp.issparse(X):
        X = X.tocsr()

    with pysam.AlignmentFile(tmp_unsorted, "wb", header=header) as bam_out:
        for cell_idx, barcode in enumerate(adata_sub.obs_names):
            if sp.issparse(X):
                row = X.getrow(cell_idx).toarray().flatten()
            else:
                row = X[cell_idx]

            for gene_idx, count in enumerate(row):
                if count <= 0:
                    continue
                gene_name = adata_sub.var_names[gene_idx]
                if gene_name not in gene_coords:
                    skipped_genes += 1
                    continue
                chrom, tss, strand = gene_coords[gene_name]
                if chrom not in chrom_sizes:
                    skipped_genes += 1
                    continue

                # Clamp to chromosome bounds
                pos = min(max(tss, 0), chrom_sizes[chrom] - read_length)

                for _ in range(int(count)):
                    a = pysam.AlignedSegment(header)
                    a.query_name = f"{barcode}_{gene_name}_{written}"
                    a.query_sequence = "N" * read_length
                    a.flag = _FLAG_REVERSE if strand == "-" else 0
                    a.reference_id = header.get_tid(chrom)
                    a.reference_start = pos
                    a.mapping_quality = 255
                    a.cigar = [(0, read_length)]  # M
                    a.query_qualities = pysam.qualitystring_to_array("I" * read_length)
                    bam_out.write(a)
                    written += 1

    print(f"[counts_to_bam] Reads written: {written}  |  genes skipped: {skipped_genes}")

    # Sort and index
    print(f"[counts_to_bam] Sorting → {output_bam}")
    pysam.sort("-o", output_bam, tmp_unsorted)
    pysam.index(output_bam)
    os.remove(tmp_unsorted)
    print(f"[counts_to_bam] Done: {output_bam}")


if __name__ == "__main__":
    try:
        _anndata_path = snakemake.input.anndata  # type: ignore[name-defined]
        _barcodes_path = snakemake.input.barcodes  # type: ignore[name-defined]
        _output_bam = snakemake.output.bam  # type: ignore[name-defined]
        _chrom_sizes = snakemake.params.chrom_sizes  # type: ignore[name-defined]
        _gtf = snakemake.params.get("gtf", "")  # type: ignore[name-defined]
        _read_length = int(snakemake.params.get("read_length", READ_LENGTH_DEFAULT))  # type: ignore[name-defined]
    except NameError:
        if len(sys.argv) < 5:
            print(
                f"Usage: {sys.argv[0]} "
                "<adata.h5ad> <barcodes.txt> <out.bam> <chrom.sizes> [gtf] [read_length]"
            )
            sys.exit(1)
        _anndata_path = sys.argv[1]
        _barcodes_path = sys.argv[2]
        _output_bam = sys.argv[3]
        _chrom_sizes = sys.argv[4]
        _gtf = sys.argv[5] if len(sys.argv) > 5 else ""
        _read_length = int(sys.argv[6]) if len(sys.argv) > 6 else READ_LENGTH_DEFAULT

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(_output_bam)), exist_ok=True)
    main(_anndata_path, _barcodes_path, _output_bam, _chrom_sizes, _gtf, _read_length)
