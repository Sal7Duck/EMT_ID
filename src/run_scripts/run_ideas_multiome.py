import sys
from pathlib import Path

import anndata as ad

# Make the local IDEAS source available
sys.path.insert(0, "/home/csb/sampoorn/IDEAS")
from computeID.functions_IDEAS import compute_ID


# Input/output locations
INPUT = Path("/home/csb/sampoorn/emt_scc/multiome_scc/Multiome_scc_IDEAS_input.h5ad")

OUTPUT_CSV = INPUT.parent / "IDEAS_local_2NN_multiome_SCC_600.csv"
OUTPUT_PKL = INPUT.parent / "IDEAS_local_2NN_multiome_SCC_600.pkl"
OUTPUT_DONE = INPUT.parent / "IDEAS_multiome_SCC_600_FINISHED.txt"


def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Input AnnData file was not found: {INPUT}")

    print("Loading AnnData...", flush=True)
    adata = ad.read_h5ad(INPUT)

    print(adata, flush=True)
    print("Shape:", adata.shape, flush=True)
    print("X type:", type(adata.X), flush=True)
    print("X dtype:", adata.X.dtype, flush=True)

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
        id_score=False,
    )

    print("IDEAS calculation finished.", flush=True)
    print(id_result, flush=True)

    id_result.to_csv(OUTPUT_CSV)
    id_result.to_pickle(OUTPUT_PKL)

    OUTPUT_DONE.write_text(
        "IDEAS local_2nn calculation completed successfully "
        "with n_neighbors=600.\n"
    )

    print(f"CSV saved: {OUTPUT_CSV}", flush=True)
    print(f"PKL saved: {OUTPUT_PKL}", flush=True)
    print(f"Completion marker saved: {OUTPUT_DONE}", flush=True)
    print("DONE.", flush=True)


if __name__ == "__main__":
    main()
