# =============================================================
# IDEAS local intrinsic dimension — GSE115943 mouse serum
# Faithful replication of Staiano et al. 2025 (maddalenastn/IDEAS)
# =============================================================
import os, re, sys, time
import numpy as np, pandas as pd, scanpy as sc
from functions_IDEAS import compute_ID

TEST = '--test' in sys.argv
KS   = [100] if TEST else [100, 200]
HOME = os.path.expanduser("~/sampoorn")
OUT  = f"{HOME}/out"; os.makedirs(OUT, exist_ok=True)

t0 = time.time()
def log(m): print(f"[{(time.time()-t0)/60:6.1f}m] {m}", flush=True)
log(f"MODE = {'TEST' if TEST else 'FULL'} | k = {KS}")

# --- 1. load -------------------------------------------------
adata = sc.read_h5ad(f"{HOME}/GSE115943_serum_counts_qc.h5ad")
adata.raw = None
log(f"loaded {adata.shape}")

# --- 2. OUR QC (not in the paper) — delete these 2 lines for a
#        literal replication with no cell filtering at all
adata = adata[(adata.obs.pct_counts_mt < 10) & (adata.obs.n_genes > 500)].copy()
log(f"after QC {adata.shape}")

# --- 3. raw counts into .X -----------------------------------
#     compute_ID normalises internally; our .X is log1p(CP10K),
#     so we must hand it the integer counts instead.
adata.X = adata.layers['counts'].copy()
del adata.layers['counts']

# --- 4. PAPER §4.1.1 gene filter -----------------------------
#     "excluded genes expressed in fewer than 1% of cells and
#      retained only protein-coding genes"  -> they got ~12,000
det  = np.asarray((adata.X > 0).sum(0)).ravel() >= 0.01 * adata.n_obs
drop = re.compile(r'(^Gm\d+$)|(Rik\d*$)|(-ps\d*$)|(^A[A-Z]\d{6}$)')
adata = adata[:, det & np.array([drop.search(g) is None
                                 for g in adata.var_names])].copy()
log(f"after gene filter {adata.shape} | dense ~{adata.n_obs*adata.n_vars*4/1e9:.1f} GB")

# --- 5. roots ------------------------------------------------
if TEST:
    roots = list(np.random.default_rng(0).choice(adata.obs_names, 50, replace=False))
    log(f"TEST: {len(roots)} roots")
else:
    roots = None            # default = every cell, as in the paper's notebook

# --- 6. the IDEAS call, once per k ---------------------------
for k in KS:
    log(f"=== n_neighbors = {k} ===")
    r = compute_ID(adata, method='local_2nn', n_neighbors=k, roots=roots,
                   normalization=True, norm_sum=1e6)
    idx  = roots if roots is not None else list(adata.obs_names)
    vals = np.asarray(r).ravel()[:len(idx)]
    tag  = f"test_k{k}" if TEST else f"k{k}"
    pd.Series(vals, index=idx, name=f'local_ID_k{k}').to_csv(f"{OUT}/local_ID_{tag}.csv")
    log(f"k={k} saved | ID {vals.min():.1f}-{vals.max():.1f} median {np.median(vals):.1f}")

# --- 7. metadata for plotting per timepoint ------------------
if not TEST:
    adata.obs.to_csv(f"{OUT}/obs.csv")
log("DONE")
