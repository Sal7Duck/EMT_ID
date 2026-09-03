import anndata as ad

INPUT = "/home/csb/sampoorn/emt_scc/12_scc/emt_scc_ideas_input.h5ad"

adata = ad.read_h5ad(INPUT)

print(adata)
print("Shape:", adata.shape)

import pandas as pd

annotation_df = pd.read_csv(
    "EMT_celltype_annotations.csv",
    index_col="cell_id"
)

# Make sure both use strings for cell IDs
annotation_df.index = annotation_df.index.astype(str)
adata.obs_names = adata.obs_names.astype(str)

# Merge annotation by cell ID
adata.obs["celltype_annotation"] = (
    annotation_df["celltype_annotation"]
    .reindex(adata.obs_names)
)

print(adata.obs["celltype_annotation"].value_counts(dropna=False))

print("Total cells:", adata.n_obs)
print(
    "Annotated cells:",
    adata.obs["celltype_annotation"].notna().sum()
)
print(
    "Missing annotations:",
    adata.obs["celltype_annotation"].isna().sum()
)

print("\nFirst 10 cells:")
print(adata.obs[["celltype_annotation"]].head(10))

import sys
import matplotlib.pyplot as plt

# IDEAS source
sys.path.insert(0, "/home/csb/sampoorn/IDEAS")

from computeID.functions_IDEAS import compute_ID


# ---------------------------------------------------------
# Balanced global 2NN ID
# ---------------------------------------------------------
# Smallest population = Hybrid = 3608
# 80% of 3608 = 2886
# Therefore every state uses 2886 cells per replicate

id_reps = compute_ID(
    adata,
    group="celltype_annotation",
    method="2nn",
    sample_size=2886,
    n_samples=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False
)

print(id_reps)
print("Shape:", id_reps.shape)
