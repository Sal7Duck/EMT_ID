import sys

# Make the local IDEAS source available
sys.path.insert(0, "/home/csb/sampoorn/IDEAS")

import anndata as ad
from computeID.functions_IDEAS import compute_ID


# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT = "/home/csb/sampoorn/emt_scc/12_scc/emt_scc_ideas_input.h5ad"

OUTPUT_CSV = "/home/csb/sampoorn/emt_scc/12_scc/IDEAS_local_2NN_SCC_600.csv"
OUTPUT_PKL = "/home/csb/sampoorn/emt_scc/12_scc/IDEAS_local_2NN_SCC_600.pkl"
OUTPUT_DONE = "/home/csb/sampoorn/emt_scc/12_scc/IDEAS_SCC_600_FINISHED.txt"


# ============================================================
# LOAD DATA
# ============================================================

print("Loading AnnData...", flush=True)

adata = ad.read_h5ad(INPUT)

print(adata, flush=True)
print("Shape:", adata.shape, flush=True)
print("X type:", type(adata.X), flush=True)
print("X dtype:", adata.X.dtype, flush=True)


# ============================================================
# RUN IDEAS
# ============================================================

print("Starting IDEAS local_2nn...", flush=True)
print("n_neighbors = 600", flush=True)

id_result = compute_ID(
    adata,
    method="local_2nn",
    roots=None,
    n_neighbors=600,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False
)


# ============================================================
# SAVE RESULTS
# ============================================================

print("IDEAS calculation finished.", flush=True)
print(id_result, flush=True)

id_result.to_csv(OUTPUT_CSV)
id_result.to_pickle(OUTPUT_PKL)

print("CSV saved:", OUTPUT_CSV, flush=True)
print("PKL saved:", OUTPUT_PKL, flush=True)


# ============================================================
# COMPLETION MARKER
# ============================================================

with open(OUTPUT_DONE, "w") as f:
    f.write(
        "IDEAS local_2nn calculation completed successfully "
        "with n_neighbors=600.\n"
    )

print("DONE.", flush=True)
