import anndata as ad
import numpy as np
import scipy.sparse as sp

# ------------------------------------------------------------
# Load the IDEAS input
# ------------------------------------------------------------

input_file = "ideas_input_GSE154989_core_timecourse.h5ad"

adata = ad.read_h5ad(input_file)

print(adata)

# ------------------------------------------------------------
# BASIC
# ------------------------------------------------------------

print("\n===== BASIC =====")
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)
print("Shape:", adata.shape)

# ------------------------------------------------------------
# EXPRESSION MATRIX
# ------------------------------------------------------------

print("\n===== X =====")
print("Type:", type(adata.X))
print("dtype:", adata.X.dtype)

if sp.issparse(adata.X):
    values = adata.X.data
else:
    values = np.asarray(adata.X).ravel()
    values = values[values != 0]

if len(values) > 0:
    print("Minimum non-zero:", values.min())
    print("Maximum:", values.max())
    print("Mean non-zero:", values.mean())
    print("Median non-zero:", np.median(values))
    print("Integer-valued:", np.allclose(values, np.round(values)))

# ------------------------------------------------------------
# METADATA
# ------------------------------------------------------------

print("\n===== OBS =====")
print("obs columns:")
print(list(adata.obs.columns))

print("\n===== VAR =====")
print("var columns:")
print(list(adata.var.columns))

# ------------------------------------------------------------
# EMBEDDINGS
# ------------------------------------------------------------

print("\n===== EMBEDDINGS =====")
print("obsm keys:", list(adata.obsm.keys()))

for key in adata.obsm.keys():
    try:
        print(f"{key}: shape = {adata.obsm[key].shape}")
    except:
        pass

# ------------------------------------------------------------
# LAYERS
# ------------------------------------------------------------

print("\n===== LAYERS =====")
print("layers:", list(adata.layers.keys()))

# ------------------------------------------------------------
# UNS
# ------------------------------------------------------------

print("\n===== UNS =====")
print("uns keys:", list(adata.uns.keys()))

# ------------------------------------------------------------
# IDENTIFIERS
# ------------------------------------------------------------

print("\n===== CELL / GENE IDs =====")
print("Cell IDs unique:", adata.obs_names.is_unique)
print("Gene IDs unique:", adata.var_names.is_unique)

# ------------------------------------------------------------
# NaN / INF
# ------------------------------------------------------------

print("\n===== NaN / INF =====")
print("NaN:", np.isnan(values).any())
print("Inf:", np.isinf(values).any())

# ------------------------------------------------------------
# SPARSITY
# ------------------------------------------------------------

print("\n===== SPARSITY =====")

if sp.issparse(adata.X):

    print("Non-zero entries:", adata.X.nnz)

    total = adata.n_obs * adata.n_vars

    print("Total entries:", total)

    print(
        "Sparsity:",
        1 - adata.X.nnz / total
    )

# ------------------------------------------------------------
# MEMORY
# ------------------------------------------------------------

print("\n===== MEMORY =====")

if sp.issparse(adata.X):

    print(
        "X data memory:",
        adata.X.data.nbytes / (1024**3),
        "GiB"
    )

else:

    print(
        "X memory:",
        adata.X.nbytes / (1024**3),
        "GiB"
    )

import sys
print(sys.executable)
print(sys.version)

import sys
sys.path.insert(0, "/home/csb/sampoorn/IDEAS")

from computeID.functions_IDEAS import compute_ID

print(compute_ID)

# Inspect the distribution of expression values

import numpy as np
import scipy.sparse as sp

X = adata.X

if sp.issparse(X):
    values = X.data
else:
    values = np.asarray(X)
    values = values[values != 0]

print("Number of non-zero values:", len(values))
print("Min:", values.min())
print("Max:", values.max())
print("Mean:", values.mean())
print("Median:", np.median(values))

print("\nQuantiles:")
print(np.quantile(values, [0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]))

id_result = compute_ID(
    adata,
    method="local_2nn",
    roots=None,
    n_neighbors=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False
)
