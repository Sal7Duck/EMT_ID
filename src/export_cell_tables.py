#!/usr/bin/env python3
"""
Build one per-cell table per dataset: barcode -> ID score + entropy + annotations.

These are the tables the figure notebooks read. Every row is keyed by the cell
barcode exactly as it appears in the .h5ad object, so a table can be joined back
onto the object with a plain index merge.

Run once, from the repo root, in the `emt-id` environment:

    python src/export_cell_tables.py --data-dir /path/to/h5ads

Writes scores/<dataset>/cell_table.csv.gz.
"""

import argparse
import os

import h5py
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:                              # anndata < 0.11
    from anndata.experimental import read_elem


def read_obs(path):
    """Read only .obs. The mouse object is 8.6 GB; read_h5ad would load all of it."""
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def keep(obs, cols):
    """Subset to the columns that are actually present, in the order given."""
    present = [c for c in cols if c in obs.columns]
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        print(f"    absent, skipped: {missing}")
    return obs[present].copy()


DATASETS = {
    # name -> (h5ad filename, obs columns to keep, extra score CSVs to join)
    "mouse_reprogramming": dict(
        h5ad="GSE115943_serum_IDEAS_with_fates_CCAT.h5ad",
        cols=[
            "timepoint", "condition", "orig_barcode",
            "IDEAS_score", "CCAT_score",
            "KS_score", "GS76_score",
            "score_ks_epithelial", "score_ks_mesenchymal", "score_pluripotency",
            "w_iPSC", "w_Stromal", "w_Trophoblast", "w_Neural",
            "P_iPSC", "P_Stromal", "P_Trophoblast", "P_Neural", "P_Other",
        ],
        join=[],
    ),
    # The ssGSEA object is the same 23,959 cells as 12_emt_scc_ID_score.h5ad with
    # HPCS_ssGSEA added, so it is the one read here. `cluster` and
    # `cluster_annotation` are TWO DIFFERENT labellings of these cells — see
    # NOTES_open_questions.md. Figure 3a uses `cluster_annotation`; the HPCS
    # programme panel uses `cluster`.
    "scc_12tumour": dict(
        h5ad="12_emt_scc_ID_score_HPCS_ssGSEA.h5ad",
        cols=[
            "cluster_annotation", "celltype_annotation", "cluster",
            "sample", "assay",
            "ID_score", "CCAT_score", "transcriptome_entropy", "HPCS_ssGSEA",
        ],
        join=[("KS_ssGSEA.csv", "KS_Epi_ssGSEA", "KS_Epi_ssGSEA"),
              ("KS_ssGSEA.csv", "KS_Mes_ssGSEA", "KS_Mes_ssGSEA")],
    ),
    "scc_multiome": dict(
        h5ad="adata_with_emt_state.h5ad",
        cols=["emt_state"],
        join=[("IDEAS_local_2NN_multiome_SCC_600.csv", "ID", "ID_score")],
    ),
    "hpcs": dict(
        h5ad="GSE154989_with_IDAS_ID_k100.h5ad",
        cols=[
            "cell_type", "lineage", "mouseID", "sampleID",
            "condition", "genotype", "sex", "week", "timesimple",
            "IDEAS_local_ID_100", "IDEAS_ID_minmax_100",
        ],
        join=[],
    ),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data",
                    help="directory holding the .h5ad objects")
    ap.add_argument("--scores-dir", default="scores")
    args = ap.parse_args()

    for name, spec in DATASETS.items():
        path = os.path.join(args.data_dir, spec["h5ad"])
        if not os.path.exists(path):
            print(f"[skip] {name}: {path} not found")
            continue

        print(f"[{name}] reading obs from {spec['h5ad']}")
        tab = keep(read_obs(path), spec["cols"])

        for fname, src_col, dst_col in spec["join"]:
            csv = os.path.join(args.scores_dir, name, fname)
            extra = pd.read_csv(csv, index_col=0)
            tab[dst_col] = extra[src_col].reindex(tab.index)
            print(f"    joined {dst_col} from {fname} "
                  f"({tab[dst_col].notna().sum():,}/{len(tab):,} matched)")

        tab.index.name = "barcode"
        out = os.path.join(args.scores_dir, name, "cell_table.csv.gz")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        # 6 significant figures: keeps every score exact to well past what is
        # reported, and keeps the mouse table small enough to sit in git.
        tab.to_csv(out, compression="gzip", float_format="%.6g")
        print(f"    -> {out}  {tab.shape[0]:,} cells x {tab.shape[1]} columns")
        print(f"       {list(tab.columns)}\n")


if __name__ == "__main__":
    main()
