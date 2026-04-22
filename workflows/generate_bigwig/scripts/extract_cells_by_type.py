"""
extract_cells_by_type.py
------------------------
Extract the barcodes belonging to a single cell type from a MuData object
and write them to a plain-text file (one barcode per line).

The script handles two output modes, controlled by ``snakemake.params.modality``:
  - "rna"  → barcodes from ``mdata.mod['rna'].obs``
  - "atac" → barcodes from ``mdata.mod['atac'].obs``

Cell-type membership is determined by the annotation column in ``mdata.obs``
(auto-detected or user-supplied via config).

Usage (via Snakemake):
    snakemake.input.muon_object   — path to .h5mu
    snakemake.output.barcodes     — output barcode text file
    snakemake.params.celltype     — target cell-type label (string)
    snakemake.params.modality     — "rna" or "atac"
    snakemake.params.celltype_column — column name, or None / "" for auto-detect
"""

import sys
import muon
import pandas


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

CANDIDATE_COLUMNS = [
    "cell_type",
    "celltype",
    "annotation",
    "leiden_celltype",
    "leiden",
]


def detect_celltype_column(obs, hint=None):
    if hint:
        if hint not in obs.columns:
            raise ValueError(
                f"Specified celltype_column '{hint}' not found.\n"
                f"Available: {list(obs.columns)}"
            )
        return hint
    for col in CANDIDATE_COLUMNS:
        if col in obs.columns:
            print(f"[extract_cells] Auto-detected cell-type column: '{col}'")
            return col
    raise ValueError(
        f"Could not auto-detect a cell-type column. Tried: {CANDIDATE_COLUMNS}\n"
        f"Available: {list(obs.columns)}"
    )


def get_barcodes_for_celltype(mdata, celltype, modality, col):
    """Return the list of barcodes for *celltype* within *modality*.

    Parameters
    ----------
    mdata : muon.MuData
    celltype : str
        Target cell-type label.
    modality : str
        "rna" or "atac".
    col : str
        Name of the cell-type annotation column in ``mdata.obs``.

    Returns
    -------
    list[str]
    """
    # Cells that belong to the requested type (global index)
    global_mask = mdata.obs[col].astype(str) == celltype
    selected_global = set(mdata.obs.index[global_mask])

    if not selected_global:
        raise ValueError(
            f"No cells found for cell type '{celltype}' "
            f"in column '{col}'."
        )

    # Intersect with modality-specific barcodes
    mod_obs = mdata.mod[modality].obs
    barcodes = [bc for bc in mod_obs.index if bc in selected_global]

    if not barcodes:
        # Fallback: try stripping modality prefix/suffix added by muon
        # (e.g. "rna:ACGT-1" vs "ACGT-1")
        stripped_global = {bc.split(":")[-1] for bc in selected_global}
        barcodes = [bc for bc in mod_obs.index if bc.split(":")[-1] in stripped_global]

    return barcodes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(muon_path, output_path, celltype, modality, hint):
    print(f"[extract_cells] Reading: {muon_path}")
    mdata = muon.read_h5mu(muon_path)

    col = detect_celltype_column(mdata.obs, hint=hint)
    barcodes = get_barcodes_for_celltype(mdata, celltype, modality, col)

    print(
        f"[extract_cells] '{celltype}' / {modality}: "
        f"{len(barcodes)} barcodes → {output_path}"
    )

    with open(output_path, "w") as fh:
        for bc in barcodes:
            fh.write(bc + "\n")


if __name__ == "__main__":
    try:
        _muon_path = snakemake.input.muon_object  # type: ignore[name-defined]
        _output_path = snakemake.output.barcodes  # type: ignore[name-defined]
        _celltype = snakemake.params.celltype  # type: ignore[name-defined]
        _modality = snakemake.params.modality  # type: ignore[name-defined]
        _hint = snakemake.params.get("celltype_column", None)  # type: ignore[name-defined]
    except NameError:
        if len(sys.argv) < 5:
            print(
                f"Usage: {sys.argv[0]} "
                "<muon.h5mu> <barcodes.txt> <celltype> <rna|atac> [column]"
            )
            sys.exit(1)
        _muon_path = sys.argv[1]
        _output_path = sys.argv[2]
        _celltype = sys.argv[3]
        _modality = sys.argv[4]
        _hint = sys.argv[5] if len(sys.argv) > 5 else None

    main(_muon_path, _output_path, _celltype, _modality, _hint)
