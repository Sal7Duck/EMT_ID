import sys

sys.path.insert(0, "/home/csb/sampoorn/IDEAS")

import anndata as ad
from computeID.functions_IDEAS import compute_ID

INPUT = "/home/csb/sampoorn/mouse_reprogramming/IDEAS_input_mouse_reprogramming.h5ad"

OUTPUT_CSV = "/home/csb/sampoorn/mouse_reprogramming/IDEAS_local_2NN_mouse_reprogramming.csv"
OUTPUT_PKL = "/home/csb/sampoorn/mouse_reprogramming/IDEAS_local_2NN_mouse_reprogramming.pkl"
OUTPUT_DONE = "/home/csb/sampoorn/mouse_reprogramming/IDEAS_FINISHED.txt"

print("Loading AnnData...", flush=True)

adata = ad.read_h5ad(INPUT)

print(adata, flush=True)
print("Shape:", adata.shape, flush=True)
print("X type:", type(adata.X), flush=True)
print("X dtype:", adata.X.dtype, flush=True)

print("Starting IDEAS local_2nn...", flush=True)

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

print("IDEAS calculation finished.", flush=True)
print(id_result, flush=True)

id_result.to_csv(OUTPUT_CSV)
id_result.to_pickle(OUTPUT_PKL)

print("CSV saved:", OUTPUT_CSV, flush=True)
print("PKL saved:", OUTPUT_PKL, flush=True)

with open(OUTPUT_DONE, "w") as f:
    f.write("IDEAS local_2nn calculation completed successfully.\n")

print("DONE.", flush=True)
