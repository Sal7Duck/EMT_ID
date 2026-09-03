#!/usr/bin/env python3
"""
Cluster-level ID: one 2NN estimate per random subsample of a labelled group.

States differ in size and the 2NN estimate depends on sample size, so every state is
subsampled to the same n, repeatedly. Needs the IDEAS environment (see environment-ideas.yml).

    python src/run_cluster_id.py --object X.h5ad --group emt_state \
        --sample-size 1377 --n-samples 100 --out scores/.../clusterID.csv

`compute_ID` takes no seed argument, so numpy's global seed is set here — without it every
run gives a different draw and the deposited table cannot be reproduced.
"""

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from computeID.functions_IDEAS import compute_ID


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--object", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--sample-size", type=int, required=True)
    ap.add_argument("--n-samples", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--state", default=None,
                    help="run one group only, writing <out>.parts/<state>.csv; "
                         "re-run per state, then --merge")
    ap.add_argument("--merge", action="store_true",
                    help="combine the per-state parts into the final table")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    parts = args.out + ".parts"

    if args.merge:
        import glob
        cols = {}
        for f in sorted(glob.glob(os.path.join(parts, "*.csv"))):
            state = os.path.splitext(os.path.basename(f))[0]
            cols[state] = pd.read_csv(f, index_col=0).iloc[:, 0]
        table = pd.DataFrame(cols)
        table.index.name = "replicate"
        table.to_csv(args.out, float_format="%.6g")
        print(f"merged {len(cols)} states -> {args.out}  {table.shape}")
        print(table.describe().round(4).to_string())
        return

    np.random.seed(args.seed)

    adata = ad.read_h5ad(args.object)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)
    print(adata.obs[args.group].value_counts().to_string(), flush=True)

    if args.state:
        # One group per invocation: each call has to finish inside the shell's time limit.
        adata = adata[adata.obs[args.group] == args.state].copy()
        adata.obs[args.group] = adata.obs[args.group].astype(str)
        print(f"state {args.state!r}: {adata.n_obs:,} cells", flush=True)

    id_reps = compute_ID(
        adata,
        group=args.group,
        method="2nn",
        sample_size=args.sample_size,
        n_samples=args.n_samples,
        normalization=True,
        norm_sum=1e6,
        full_output=True,
        id_score=False,
    )

    id_reps.index.name = "replicate"
    if args.state:
        os.makedirs(parts, exist_ok=True)
        out = os.path.join(parts, f"{args.state}.csv")
    else:
        out = args.out
        os.makedirs(os.path.dirname(out), exist_ok=True)
    id_reps.to_csv(out, float_format="%.6g")
    args.out = out
    print(f"\nwritten: {args.out}  {id_reps.shape}", flush=True)
    print(id_reps.describe().round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
