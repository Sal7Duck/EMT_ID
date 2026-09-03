# Figure_scc_multiome

Everything needed to reproduce the SCC multiome EMT-state figure, and a record of how every
score and label in it was generated.

```
Figure_scc_multiome/
  README.md
  Figure_scc_multiome.ipynb          the notebook — three panels
  Input_scc_multiome_clusterID.csv   100 replicates x 5 states — the ID panel
  Input_scc_multiome_percell.csv     one row per barcode — the KS scatter
  figures/                           600 dpi PNG + vector PDF
  source_data/                       the numbers behind each panel
```

> **The ID in this figure is cluster-level.** IDEAS `method="2nn"` returns one ID per random
> subsample of a state, not one per cell. A per-cell ID also exists for these cells
> (`ID_score_percell`) and is **not** what the ID panel plots.

Open the notebook and run it top to bottom. Environment: `../environment.yml` (`emt-id`).

---

## Where the data came from

**[GSE288220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288220)** — "Single cell
multiomics unravel the transcription networks controlling the different EMT tumor states
[multiome]", mouse squamous cell carcinoma.

**This figure uses the RNA modality only.** The multiome assay measures gene expression and
chromatin accessibility in the same nuclei; everything here — the intrinsic dimension, the KS
signature scores — is computed on the RNA side. The ATAC modality is not used.

Source object: **`adata_with_emt_state.h5ad`**, 14,302 cells × 26,408 genes, raw integer
counts, mouse gene symbols, carrying the five EMT states in `obs["emt_state"]`:

| State | Cells |
|---|---|
| Epithelial | 5,827 |
| Early hybrid (EMT-I) | 2,472 |
| Early hybrid (EMT-II) | 1,830 |
| Hybrid | 1,722 |
| Late EMT | 2,451 |

---

## Why there are two input files

The two panels sit on **different units**.

The **ID panel is cluster-level**: `compute_ID` with `method="2nn"` and `group=` returns one
ID per random subsample of a state, not one per cell. Each of the 100 numbers per state is a
property of a subsample, so there is no barcode to attach it to and it cannot be merged into
`obs`. It gets its own file.

The **KS scatter is per-cell**, so it lives in the barcode-indexed table.

A per-cell ID also exists for this dataset (`local_2nn`, k = 600, in
`../scores/scc_multiome/`) and is carried in the per-cell table as `ID_score_percell` — but
it is **not** what this figure's ID panel plots.

---

## How each score was generated

### The cluster-level ID

The five states range from 1,722 to 5,827 cells. The 2NN estimate depends on sample size, so
states are compared at equal *n*: **1,377 cells** — 80 % of the smallest state — drawn at
random, **100 times per state**.

**IDEAS** — <https://github.com/maddalenastn/IDEAS>, vendored in `../src/computeID/`.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922.
> [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)

```python
id_reps = compute_ID(
    adata,                          # RNA modality, raw counts
    group="emt_state",
    method="2nn",
    sample_size=1377,
    n_samples=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,
)                                   # -> 100 x 5
```

Originally run by `Multiome_cluster_id.ipynb`, where the result was never saved. Re-run here
by `../src/run_cluster_id.py` **with `--seed 0`**, so the deposited table is reproducible:

```bash
for S in "Epithelial" "Early hybrid (EMT-I)" "Early hybrid (EMT-II)" "Hybrid" "Late EMT"; do
  python src/run_cluster_id.py \
    --object adata_with_emt_state.h5ad --group emt_state \
    --sample-size 1377 --n-samples 100 --seed 0 --state "$S" \
    --out scores/scc_multiome/IDEAS_clusterID_100runs.csv
done
python src/run_cluster_id.py --merge --object x --group emt_state --sample-size 1377 \
  --out scores/scc_multiome/IDEAS_clusterID_100runs.csv
```

(One state per call, because each takes about 100 seconds.)

Median raw ID per state, which is what the slide-deck version plots:

| State | Median raw ID |
|---|---|
| Early hybrid (EMT-I) | 31.49 |
| Epithelial | 28.70 |
| Early hybrid (EMT-II) | 27.53 |
| Late EMT | 26.56 |
| Hybrid | 26.09 |

The input file carries both: `<state>` on 0–1 and `<state>__raw` as `compute_ID` returned it.

### `KS_Epi_ssGSEA`, `KS_Mes_ssGSEA`

ssGSEA enrichment scores for the mouse KS epithelial (164 genes) and KS mesenchymal (48
genes) signatures, `gseapy.ssgsea` with `sample_norm_method="rank"`, on the same RNA matrix.
Computed inside `Multiome_cluster_id.ipynb` and never saved with the object;
`../src/compute_ks_ssgsea.py` recomputes them and writes
`../scores/scc_multiome/KS_ssGSEA.csv`. Reproduces the published r = −0.074, p = 5.21e-19,
n = 14,302.

---

## Scaling

All scores 0–1. The cluster-level ID is min–max scaled across the whole 100 × 5 table; the
per-cell scores across all 14,302 cells. Raw values are kept alongside in `_raw` / `__raw`
columns.

Min–max is linear, so Kruskal-Wallis *H*, Dunn *q* and Pearson *r* are unchanged.

---

## References

| | |
|---|---|
| Dataset | [GSE288220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288220) — RNA modality of the multiome |
| IDEAS | <https://github.com/maddalenastn/IDEAS> · [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
| gseapy | <https://github.com/zqfang/GSEApy> |
