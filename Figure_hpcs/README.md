# Figure_hpcs

Everything needed to reproduce the HPCS figure, and a record of how every score and label in
it was generated.

```
Figure_hpcs/
  README.md
  Figure_hpcs.ipynb              the notebook — two panels
  Input_hpcs_percell.csv         3,891 barcodes — lung tumour cells
  figures/                       600 dpi PNG + vector PDF
  source_data/                   the numbers behind each panel
```

Open the notebook and run it top to bottom. Environment: `../environment.yml` (`emt-id`).

---

## Where the data came from

**[GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989)** — "Emergence of
a high-plasticity cell state during lung cancer evolution, single-cell RNAseq timecourse".
KP/KPT mouse lung adenocarcinoma. Source object `GSE154989_with_IDAS_ID_k100.h5ad` — 3,891
cells, 12 annotated cell types, 39 mice.

The **HPCS programme scored across the SCC EMT states** is not here. That panel is measured
on SCC cells and grouped by EMT state, so it lives with the rest of that cohort in
`../Figure_scc_12tumour/` (panels 4–5).

---

## How each score was generated

### `ID_score` — intrinsic dimension

Per-cell IDEAS, `method="local_2nn"`, **k = 100**.

**IDEAS** — <https://github.com/maddalenastn/IDEAS>, vendored in `../src/computeID/`.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922.
> [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)

```python
id_result = compute_ID(
    adata,                      # raw counts
    method="local_2nn",
    roots=None,
    n_neighbors=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,
)
```

Run by `../src/run_scripts/run_hpcs_ideas.py`. Raw output:
`../scores/hpcs/IDEAS_local_2NN_GSE154989_k100.csv`.

**k is 100 here and 600 everywhere else**, because this dataset has only 3,891 cells. Local
ID rises with k, so ID is compared *within* a dataset and never across datasets. This is an
open item.

---

## Statistics

Kruskal-Wallis across groups, then Dunn's post-hoc with pooled ranks, tie correction and
Benjamini-Hochberg FDR. Effect sizes are median differences and rank-biserial correlation.

**Panel 2 result:** HPCS has a higher median ID than **10 of the other 11 cell types** at
FDR < 0.05. The exception is **EMT / fibroblast-like (Vim/Axl/CD109)**, q = 0.095 — the one
state HPCS is not distinguishable from. Worth stating plainly rather than rounding to
"higher than all".

---

## One thing to know about the numbers

**The unit of replication is the mouse, not the cell.** All tests here treat cells as
independent; cells from one animal are not. The per-cell table carries `mouseID` (39 animals)
so a mouse-level analysis can be run. At mouse level HPCS beats 10/11 cell types and
EMT/fibroblast-like is again the exception.

---

## References

| | |
|---|---|
| Lung tumour time course | [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989) |
| IDEAS | <https://github.com/maddalenastn/IDEAS> · [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
