# Figure_scc_12tumour

Everything needed to reproduce the SCC EMT-state figure, and a record of how every score and
label in it was generated.

```
Figure_scc_12tumour/
  README.md
  Figure_scc_12tumour.ipynb          the notebook — five panels
  Input_scc_12tumour_clusterID.csv   100 replicates x 4 states — the ID panel
  Input_scc_12tumour_percell.csv     one row per barcode — the KS scatter
  figures/                           600 dpi PNG + vector PDF
  source_data/                       the numbers behind each panel
```

> **The ID in this figure is cluster-level.** IDEAS `method="2nn"` returns one ID per random
> subsample of a state, not one per cell. A per-cell ID also exists for these cells
> (`ID_score_percell`) and is **not** what the ID panel plots.

Open the notebook and run it top to bottom. Environment: `../environment.yml` (`emt-id`).

---

## Dataset

**[GSE288137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288137)** — mouse
squamous cell carcinoma, 12 tumours, 23,959 cells (14,059 scRNA + 9,900 multiome), four EMT
states.

---

## Why there are two input files

The two panels sit on **different units**, and this is the thing to understand about the
figure.

The **ID panel is cluster-level**. IDEAS was run with `method="2nn"` and `group=`, which
returns one ID per random subsample of a state — not one per cell. Each of the 100 numbers
per state is a property of a subsample, so there is no barcode to attach it to and it cannot
be merged into `obs`. It gets its own file.

The **KS scatter is per-cell**, so it lives in the barcode-indexed table.

A per-cell ID does exist for this dataset (`local_2nn`, k = 600, in
`../scores/scc_12tumour/`) and is carried in the per-cell table as `ID_score_percell` — but
it is **not** what this figure's ID panel plots.

---

## How each score was generated

### The cluster-level ID

The four states differ in size from 3,608 to 9,985 cells. The 2NN estimate depends on sample
size, so states are compared at equal *n*: **2,886 cells** — 80 % of the smallest state —
drawn at random, **100 times per state**.

**IDEAS** — <https://github.com/maddalenastn/IDEAS>, vendored in `../src/computeID/`.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922.
> [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)

```python
id_reps = compute_ID(
    adata,                          # raw counts
    group="celltype_annotation",
    method="2nn",
    sample_size=2886,
    n_samples=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,
)                                   # -> 100 x 4
```

Run by `local_id.ipynb` in the `ideas` environment; about 15 minutes.

### `KS_Epi_ssGSEA`, `KS_Mes_ssGSEA`

ssGSEA enrichment scores for the mouse KS epithelial (164 genes) and KS mesenchymal (48
genes) signatures, `gseapy.ssgsea` with `sample_norm_method="rank"`. Computed inside
`id_entropy_12_emt_scc.ipynb` and never saved with the deposited object;
`../src/compute_ks_ssgsea.py` recomputes them and writes
`../scores/scc_12tumour/KS_ssGSEA.csv`. Reproduces the published r = −0.627.

### `HPCS_ssGSEA` — the HPCS programme

ssGSEA enrichment score for the HPCS signature, carried in `adata.obs["HPCS_ssGSEA"]` of the
SCC object. The signature was defined in the KP lung tumours (`../Figure_hpcs/`); panels 4–5
ask whether it also tracks EMT state here.

On `cluster_annotation` the omnibus test gives **H = 1308.36** and all six pairwise
comparisons are significant, Early-Hybrid highest and Late-EMT lowest. (The originally
published number, H = 1458.58, used the `cluster` labelling, which this folder no longer
carries.)

### `signalling_entropy`

CCAT from `SCENT::CompCCAT`. Full provenance in
`../Figure_mouse_reprogramming/README.md`.

> Teschendorff AE, Enver T. *Nat Commun* 2017;8:15599.
> [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)
> Code used: <https://github.com/Zaoqu-Liu/SCENT> (Rcpp fork of <https://github.com/aet21/SCENT>)

---

## State labels

Both panels use `cluster_annotation`: Epithelial 4,146 · Early-Hybrid 9,985 · Hybrid 3,608 ·
Late-EMT 6,220. This is the labelling `compute_ID` was grouped on.

---

## Scaling

All scores 0–1. The cluster-level ID file is IDEAS' own min–max output; the per-cell scores
are min–max scaled over all 23,959 cells, with the raw value kept in a `_raw` column.

Min–max is linear, so Kruskal-Wallis *H*, Dunn *q* and Pearson *r* are unchanged.

---

## Two things to know about the numbers

**The deposited replicate table is a different random draw from the one in the slide deck.**
`compute_ID` takes no seed, so each run of the 100 subsamples gives different numbers. The
draw behind the slides gives Kruskal-Wallis H = 313.61; the deposited draw gives H = 315.17.
Every conclusion is the same, and Epithelial vs Early-Hybrid is non-significant in both — but
they are not the same 400 numbers. **Fix before submission:** set `np.random.seed()` before
the call and deposit that draw.

**The p-values measure estimator stability, not biological replication.** The 100 replicates
resample the same 23,959 cells, so raising `n_samples` lowers every p-value. Rank-biserial is
1.000 on four of the six pairs — complete separation of the resampling distributions. The
unit of biological replication here is the tumour (`sample`, n = 12), not the subsample.

---

## References

| | |
|---|---|
| Dataset | [GSE288137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288137) |
| IDEAS | <https://github.com/maddalenastn/IDEAS> · [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
| SCENT / CCAT | <https://github.com/Zaoqu-Liu/SCENT> · <https://github.com/aet21/SCENT> · [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599) |
| gseapy | <https://github.com/zqfang/GSEApy> |
