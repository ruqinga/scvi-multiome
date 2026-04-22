"""
fragments_to_bam.py
-------------------
Convert a subset of ATAC fragment records (BED-like file) into a
coordinate-sorted, indexed BAM file for use with ``bamCoverage``.

Each line of a Signac/Cell Ranger ATAC fragment file has the format:
    chrom  start  end  barcode  count

For every fragment record whose barcode is in the supplied list the script
emits two synthetic paired-end reads anchored at the fragment ends:
  - Read 1: forward strand, starts at ``start``
  - Read 2: reverse strand, ends at ``end``

Both reads receive the same mapping quality (255) and a simple CIGAR
reflecting the fragment half-length (capped at ``read_length``).

Usage (via Snakemake)
~~~~~~~~~~~~~~~~~~~~~
    snakemake.input.fragments    — path to fragment BED (may be .gz)
    snakemake.input.barcodes     — text file, one barcode per line
    snakemake.output.bam         — output sorted + indexed BAM
    snakemake.params.chrom_sizes — path to .chrom.sizes file
    snakemake.params.read_length — synthetic read length (default 50)
"""

import sys
import os
import gzip

import pysam

# SAM flag constants (pysam >= 0.16 exposes these, but define them explicitly
# for clarity and compatibility with older versions)
_FLAG_PAIRED      = 0x1   # read is paired
_FLAG_PROPER_PAIR = 0x2   # read mapped in proper pair
_FLAG_REVERSE     = 0x10  # read mapped on reverse strand
_FLAG_READ1       = 0x40  # this is read 1
_FLAG_READ2       = 0x80  # this is read 2

# Read 1 forward:  paired + proper_pair + read1       = 0x1 | 0x2 | 0x40 = 67
_FLAG_R1_FWD = _FLAG_PAIRED | _FLAG_PROPER_PAIR | _FLAG_READ1
# Read 2 reverse:  paired + proper_pair + reverse + read2 = 0x1 | 0x2 | 0x10 | 0x80 = 131
_FLAG_R2_REV = _FLAG_PAIRED | _FLAG_PROPER_PAIR | _FLAG_REVERSE | _FLAG_READ2


READ_LENGTH_DEFAULT = 50


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def build_bam_header(chrom_sizes):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": length} for chrom, length in chrom_sizes.items()],
    }
    return pysam.AlignmentHeader.from_dict(header)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_fragment_file(bed_path, barcodes, barcodes_stripped, chrom_sizes, header,
                          bam_out, read_length, read_counter):
    """Write synthetic paired-end reads from one fragment file into *bam_out*.

    Returns (written, skipped) counts.
    """
    written = 0
    skipped = 0
    with open_maybe_gzip(bed_path) as bed_in:
        for line in bed_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom, start_s, end_s, barcode = parts[:4]

            if barcode not in barcodes and barcode not in barcodes_stripped:
                continue

            if chrom not in chrom_sizes:
                skipped += 1
                continue

            start = int(start_s)
            end = int(end_s)
            frag_len = end - start

            if frag_len <= 0:
                skipped += 1
                continue

            rl = min(read_length, frag_len)
            idx = read_counter[0]

            # Read 1 — forward, starts at fragment start
            r1 = pysam.AlignedSegment(header)
            r1.query_name = f"{barcode}_{chrom}_{start}_{end}_{idx}"
            r1.query_sequence = "N" * rl
            r1.flag = _FLAG_R1_FWD   # paired, proper pair, read1, forward
            r1.reference_id = header.get_tid(chrom)
            r1.reference_start = start
            r1.mapping_quality = 255
            r1.cigar = [(0, rl)]
            r1.query_qualities = pysam.qualitystring_to_array("I" * rl)
            r1.next_reference_id = r1.reference_id
            r1.next_reference_start = max(0, end - rl)
            r1.template_length = frag_len
            bam_out.write(r1)

            # Read 2 — reverse, ends at fragment end
            r2 = pysam.AlignedSegment(header)
            r2.query_name = r1.query_name
            r2.query_sequence = "N" * rl
            r2.flag = _FLAG_R2_REV   # paired, proper pair, read2, reverse
            r2.reference_id = r1.reference_id
            r2.reference_start = max(0, end - rl)
            r2.mapping_quality = 255
            r2.cigar = [(0, rl)]
            r2.query_qualities = pysam.qualitystring_to_array("I" * rl)
            r2.next_reference_id = r1.reference_id
            r2.next_reference_start = start
            r2.template_length = -frag_len
            bam_out.write(r2)

            written += 2
            read_counter[0] += 2

    return written, skipped


def main(fragments_paths, barcodes_path, output_bam, chrom_sizes_path, read_length):
    """
    Parameters
    ----------
    fragments_paths : str or list[str]
        One or more fragment BED (optionally .gz) file paths.
    """
    if isinstance(fragments_paths, str):
        fragments_paths = [fragments_paths]

    print(f"[fragments_to_bam] Reading barcodes: {barcodes_path}")
    with open(barcodes_path) as fh:
        barcodes = {line.strip() for line in fh if line.strip()}
    barcodes_stripped = {bc.split(":")[-1] for bc in barcodes}
    print(f"[fragments_to_bam] {len(barcodes)} barcodes selected")

    chrom_sizes = parse_chrom_sizes(chrom_sizes_path)
    header = build_bam_header(chrom_sizes)

    tmp_unsorted = output_bam + ".unsorted.bam"
    total_written = 0
    total_skipped = 0
    read_counter = [0]  # mutable counter shared across calls

    with pysam.AlignmentFile(tmp_unsorted, "wb", header=header) as bam_out:
        for frag_path in fragments_paths:
            print(f"[fragments_to_bam] Processing: {frag_path}")
            w, s = process_fragment_file(
                frag_path, barcodes, barcodes_stripped,
                chrom_sizes, header, bam_out, read_length, read_counter
            )
            total_written += w
            total_skipped += s

    print(f"[fragments_to_bam] Reads written: {total_written}  |  records skipped: {total_skipped}")

    print(f"[fragments_to_bam] Sorting → {output_bam}")
    pysam.sort("-o", output_bam, tmp_unsorted)
    pysam.index(output_bam)
    os.remove(tmp_unsorted)
    print(f"[fragments_to_bam] Done: {output_bam}")


if __name__ == "__main__":
    try:
        # Snakemake may pass a list or a single string for input.fragments
        _raw_frags = snakemake.input.fragments  # type: ignore[name-defined]
        _fragments_paths = _raw_frags if isinstance(_raw_frags, list) else [_raw_frags]
        _barcodes_path = snakemake.input.barcodes  # type: ignore[name-defined]
        _output_bam = snakemake.output.bam  # type: ignore[name-defined]
        _chrom_sizes = snakemake.params.chrom_sizes  # type: ignore[name-defined]
        _read_length = int(snakemake.params.get("read_length", READ_LENGTH_DEFAULT))  # type: ignore[name-defined]
    except NameError:
        if len(sys.argv) < 5:
            print(
                f"Usage: {sys.argv[0]} "
                "<fragments.bed[.gz]> <barcodes.txt> <out.bam> <chrom.sizes> [read_length]"
            )
            sys.exit(1)
        _fragments_paths = [sys.argv[1]]
        _barcodes_path = sys.argv[2]
        _output_bam = sys.argv[3]
        _chrom_sizes = sys.argv[4]
        _read_length = int(sys.argv[5]) if len(sys.argv) > 5 else READ_LENGTH_DEFAULT

    os.makedirs(os.path.dirname(os.path.abspath(_output_bam)), exist_ok=True)
    main(_fragments_paths, _barcodes_path, _output_bam, _chrom_sizes, _read_length)
