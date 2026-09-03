import sys

sys.path.insert(0, "/home/csb/sampoorn/IDEAS")

import anndata as ad

from computeID.functions_IDEAS import compute_ID

import anndata as ad
from computeID.functions_IDEAS import compute_ID

INPUT = "/home/csb/sampoorn/emt_scc/GSE288137_SCO_IDEAS_092.h5ad"

OUTPUT_PKL = "/home/csb/sampoorn/emt_scc/IDEAS_local_2NN_result.pkl"
OUTPUT_CSV = "/home/csb/sampoorn/emt_scc/IDEAS_local_2NN_result.csv"

print("Loading AnnData...", flush=True)

adata = ad.read_h5ad(INPUT)

print(adata, flush=True)
print("Starting IDEAS local_2nn...", flush=True)

id_result = compute_ID(
    adata,
    method="local_2nn",
    roots=None,
    n_neighbors=0,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False
)

print("IDEAS calculation finished.", flush=True)
print(id_result, flush=True)

# Save the exact returned DataFrame
id_result.to_pickle(OUTPUT_PKL)
id_result.to_csv(OUTPUT_CSV)

print("Results saved.", flush=True)
print("Pickle:", OUTPUT_PKL, flush=True)
print("CSV:", OUTPUT_CSV, flush=True)
print("IDEAS calculation finished.", flush=True)
print(id_result, flush=True)

# Save the exact returned result
id_result.to_pickle(OUTPUT_PKL)
id_result.to_csv(OUTPUT_CSV)

print("Results saved.", flush=True)
print(f"Pickle: {OUTPUT_PKL}", flush=True)
print(f"CSV: {OUTPUT_CSV}", flush=True)

# Completion marker
with open("/home/csb/sampoorn/emt_scc/IDEAS_FINISHED.txt", "w") as f:
    f.write("IDEAS local_2nn calculation completed successfully.\n")
