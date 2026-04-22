"""
get_celltypes.py
----------------
Inspect a MuData (.h5mu) object and:
  1. Auto-detect (or validate) the cell-type annotation column.
  2. Write the list of unique cell-type labels to a plain-text file
     (one label per line), so that Snakemake can use a checkpoint to
     dynamically expand output targets.

Usage (via Snakemake):
    snakemake.input.muon_object  — path to .h5mu
    snakemake.output.celltypes   — output text file
    snakemake.params.celltype_column — column name, or None / "" for auto-detect
"""

import sys
import muon


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
    """Return the first matching cell-type column name.

    Parameters
    ----------
    obs : pandas.DataFrame
        ``mdata.obs`` from a MuData object.
    hint : str or None
        User-supplied column name.  If *None* or empty, auto-detection is used.

    Returns
    -------
    str
        Column name.

    Raises
    ------
    ValueError
        When the hinted column is absent, or no candidate column is found.
    """
    if hint:
        if hint not in obs.columns:
            raise ValueError(
                f"Specified celltype_column '{hint}' not found in mdata.obs.\n"
                f"Available columns: {list(obs.columns)}"
            )
        return hint

    for col in CANDIDATE_COLUMNS:
        if col in obs.columns:
            print(f"[get_celltypes] Auto-detected cell-type column: '{col}'")
            return col

    raise ValueError(
        "Could not auto-detect a cell-type column in mdata.obs.\n"
        f"Tried: {CANDIDATE_COLUMNS}\n"
        f"Available columns: {list(obs.columns)}\n"
        "Set 'celltype_column' explicitly in config/config.yaml."
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(muon_path, output_path, celltype_column_hint):
    print(f"[get_celltypes] Reading: {muon_path}")
    mdata = muon.read_h5mu(muon_path)

    col = detect_celltype_column(mdata.obs, hint=celltype_column_hint)

    unique_labels = sorted(mdata.obs[col].dropna().astype(str).unique().tolist())

    if not unique_labels:
        raise ValueError(f"Column '{col}' contains no non-null values.")

    print(f"[get_celltypes] Found {len(unique_labels)} cell type(s): {unique_labels}")

    with open(output_path, "w") as fh:
        for label in unique_labels:
            fh.write(label + "\n")

    print(f"[get_celltypes] Written to: {output_path}")


if __name__ == "__main__":
    # Snakemake injects the `snakemake` object into the global scope when
    # the script is executed via `script:` directive.
    try:
        _muon_path = snakemake.input.muon_object  # type: ignore[name-defined]
        _output_path = snakemake.output.celltypes  # type: ignore[name-defined]
        _hint = snakemake.params.get("celltype_column", None)  # type: ignore[name-defined]
    except NameError:
        # Allow standalone execution for testing
        if len(sys.argv) < 3:
            print(f"Usage: {sys.argv[0]} <muon_object.h5mu> <output_celltypes.txt> [column_hint]")
            sys.exit(1)
        _muon_path = sys.argv[1]
        _output_path = sys.argv[2]
        _hint = sys.argv[3] if len(sys.argv) > 3 else None

    main(_muon_path, _output_path, _hint)
