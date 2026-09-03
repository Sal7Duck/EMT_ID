# Figure_mouse_reprogramming

Everything needed to reproduce the mouse reprogramming figure, and a record of how every
score and every label in it was generated.

```
Figure_mouse_reprogramming/
  README.md                          this file
  Figure_mouse_reprogramming.ipynb   the notebook — all thirteen panels
  Input_mouse_reprogramming.csv      the single input: one row per cell barcode
  figures/                           every panel, 600 dpi PNG + vector PDF
  source_data/                       the numbers behind every panel
```

**To reproduce:** open the notebook and run it top to bottom. It reads one CSV and writes
`figures/` and `source_data/`. No `.h5ad` is opened, nothing is downloaded, and nothing is
recomputed — the ID and entropy runs happen upstream (see below).

Environment: `../environment.yml` (`emt-id`) — python, pandas, numpy, scipy, statsmodels,
matplotlib, seaborn.

---

## Dataset

**[GSE115943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115943)** — mouse
embryonic fibroblasts reprogrammed by OSKM induction, serum arm. 181,124 cells, 39 time
points, day 0 to day 18.

> Schiebinger G, Shu J, Tabaka M, Cleary B, Subramanian V, Solomon A, Gould J, Liu S, Lin S,
> Berube P, Lee L, Chen J, Brumbaugh J, Rigollet P, Hochedlinger K, Jaenisch R, Regev A,
> Lander ES. Optimal-transport analysis of single-cell gene expression identifies
> developmental trajectories in reprogramming. *Cell* 2019;176(4):928–943.e22.
> [doi:10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006) · PMID 30712874

Counts were downloaded from GEO, quality-controlled, and 3,000 highly variable genes
selected. That object, holding **raw counts**, is the input to everything below.

---

## How each score was generated

### `ID_score` — intrinsic dimension

**IDEAS** — <https://github.com/maddalenastn/IDEAS>, vendored in `../src/computeID/`.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922, posted 24 July 2025.
> [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)
 Estimator `local_2nn`: a k-nearest-neighbour graph is built over all
cells, and for each cell the two-nearest-neighbour estimator of Facco et al. is applied
inside that cell's own neighbourhood of **k = 600** cells. One ID value per cell.

```python
from computeID.functions_IDEAS import compute_ID

id_result = compute_ID(
    adata,                  # raw counts — IDEAS normalises to 1e6 internally
    method="local_2nn",
    roots=None,
    n_neighbors=600,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,         # return the raw ID, not IDEAS' own 0-1 scaling
)
```

Run by `../src/run_scripts/run_mouse_ideas.py`, in the separate `ideas` environment
(`../environment-ideas.yml`). Slow — 181,124 cells at k = 600. Raw output:
`../scores/mouse_reprogramming/IDEAS_local_2NN_mouse_reprogramming.csv`.

Because `id_score=False`, the values here are the raw estimator output, 7.82 to 51.31. All
0–1 scaling is done afterwards, in this folder.

### `signalling_entropy` — CCAT

**CCAT** — the fast approximation to signalling entropy — computed with `SCENT::CompCCAT`
via `Rscript`, on the log-normalised matrix against a protein–protein interaction network.
Wrapper in `../src/scent/`.

> Method: Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
> differentiation potency from a cell's transcriptome. *Nat Commun* 2017;8:15599.
> [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)

Implementation used: <https://github.com/Zaoqu-Liu/SCENT> — an Rcpp-optimised fork of the
authors' original package, <https://github.com/aet21/SCENT>. Provenance recorded in
`adata.uns["SCENT_CCAT"]`:

| | |
|---|---|
| version / commit | 2.0.0 / `8a3a80953befa50b4320249d28e1ce0afb6c5e06` |
| method | `CompCCAT` |
| PPI network | `net13Jun12` |
| genes used | 15,225 |
| orthology | mouse symbol → human Entrez via MyGene/HomoloGene, strict 1:1 |
| input layer | `counts` |

### `KS_score`, `KS_epithelial`, `KS_mesenchymal`

The KS (Kolmogorov–Smirnov) EMT scoring scheme. `KS_epithelial` and `KS_mesenchymal` are the
two arms scored separately with `scanpy.tl.score_genes`; `KS_score` is the combined score.

### `GS76_score`

The 76-gene signature EMT score (76GS), scored with `scanpy.tl.score_genes`.

### `pluripotency`

Pluripotency signature score, `scanpy.tl.score_genes`.

### `w_*` — Waddington-OT trajectory weights

**Waddington-OT** (`wot.ot.OTModel`), ε = 0.05, λ1 = 1, λ2 = 50, 3 growth iterations,
`local_pca=30`. Transport maps were computed between every pair of adjacent time points, then
pulled back from four terminal fates defined at day 18 — **iPSC, Stromal, Trophoblast,
Neural** — each taken as the top 500 cells by its fate signature.

`w_iPSC` and friends are **ancestor trajectory weights**: within one day they sum to 1 for
each fate. They answer *which cells at this day are ancestors of that fate*, and they are the
weights to average a score against when following a trajectory back in time.

Method, parameters and gene sets in `../wot/README.md`.

### `P_*` — Waddington-OT fate probabilities

From the same transport maps, but normalised the other way: for each cell, `P_iPSC`,
`P_Stromal`, `P_Trophoblast`, `P_Neural` and `P_Other` sum to 1. They answer *where is this
cell going*.

`w` and `P` are different objects and are not interchangeable. The trajectory panels use `w`.

---

## The input file

`Input_mouse_reprogramming.csv` — 181,124 rows × 25 columns, index = cell barcode as it
appears in the annotated object.

| Column | |
|---|---|
| `day` | collection day, 0 to 18 |
| `condition` | culture condition |
| `ID_score`, `ID_score_raw` | intrinsic dimension: 0–1 scaled, and raw |
| `signalling_entropy`, `signalling_entropy_raw` | CCAT: 0–1 scaled, and raw |
| `KS_score`, `GS76_score` | + `_raw` for each |
| `KS_epithelial`, `KS_mesenchymal` | + `_raw` for each |
| `pluripotency` | + `_raw` |
| `w_iPSC`, `w_Stromal`, `w_Trophoblast`, `w_Neural` | trajectory weights |
| `P_iPSC`, `P_Stromal`, `P_Trophoblast`, `P_Neural`, `P_Other` | fate probabilities |

Derived from `GSE115943_serum_IDEAS_with_fates_CCAT.h5ad`.

The `_raw` columns are the values exactly as the tool that produced them returned them. The
unsuffixed columns are those min–max scaled to 0–1 over all 181,124 cells.

---

## Scaling

**Every score in every panel is on 0–1.** Each panel min–max scales over the cells it
actually plots, so the panel spans its own full range.

This changes nothing statistical. Min–max is a linear transform, and one-way ANOVA *F*,
Tukey-adjusted *p*, Pearson *r*, Spearman ρ and weighted means are all invariant under it.
Every number reported in the figure is identical to what the raw scales give:

| | |
|---|---|
| ID score, ANOVA across 8 time bins | F = 5114.062 |
| Signalling entropy, ANOVA across 8 time bins | F = 2520.955 |
| KS epithelial vs KS mesenchymal, days 0–7.5 | r = −0.153 |
| KS vs 76GS, days 0–7.5 | r = −0.319 |
| Signalling entropy vs ID, all cells | r = +0.277 |
| ID vs pluripotency, days ≥ 10.5 | r = +0.609 |
| KS epithelial vs KS mesenchymal, days ≥ 10.5 | r = −0.211 |
| KS vs 76GS, days ≥ 10.5 | r = −0.696 |

Axes are drawn 0–1 wherever the plotted values span that range. The two trajectory panels are
the exception: the scores are 0–1, but the weighted means occupy a narrow band inside it, so
those axes fit the data — forcing 0–1 there would flatten the panel and hide the separation
between fates, which is the point of the panel.

---

## Naming

Two terms are used throughout, on every axis and in every title:

- **ID score** — the intrinsic dimension
- **Signalling entropy** — the CCAT score

---

## Statistics

**Across ordered time bins** — one-way ANOVA as the omnibus test, Tukey HSD for the pairwise
comparisons. Time bins pool adjacent half-day points into eight bins covering days 0 to 7.5.

**Correlations** — Pearson *r* on cells with both scores present, with the OLS line through
the same cells; Spearman ρ reported alongside in `source_data/correlations.csv`. At these
sample sizes the p-value underflows the double-precision floor and is reported as
`p < 1e-300`.

**Trajectory means** — weighted by the WOT ancestor weight `w`, over cells at that day with
non-zero weight. `n_eff = (Σw)² / Σw²` is the effective sample size behind each point and is
in `source_data/wot_trajectories.csv`.

---

## Output

Every panel: `figures/panelNN_*.png` at **600 dpi** and `figures/panelNN_*.pdf` as vector,
with TrueType fonts embedded so text stays editable in Illustrator or Inkscape. Scatter
points are rasterised inside the PDF (100k+ points as vector would not open); axes, text and
lines stay vector.

| # | Panel | Cells |
|---|---|---|
| 1 | ID score across reprogramming time | days 0–7.5 |
| 2 | KS epithelial vs KS mesenchymal | days 0–7.5 |
| 3 | KS score vs 76GS score | days 0–7.5 |
| 4 | Signalling entropy vs ID score | all cells |
| 5 | ID score: pairwise Tukey HSD | days 0–7.5 |
| 6 | Signalling entropy across reprogramming time | days 0–7.5 |
| 7 | Signalling entropy: pairwise Tukey HSD | days 0–7.5 |
| 8 | ID score along the Waddington-OT trajectories | days ≥ 10.5 |
| 9 | Signalling entropy along the Waddington-OT trajectories | days ≥ 10.5 |
| 10 | ID score vs pluripotency score | days ≥ 10.5 |
| 11 | KS epithelial vs KS mesenchymal | days ≥ 10.5 |
| 12 | KS score vs 76GS score | days ≥ 10.5 |
| 13 | KS epithelial vs KS mesenchymal along the trajectories | days ≥ 10.5 |


---

## References

| | |
|---|---|
| Dataset | [GSE115943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115943) |
| Dataset paper | Schiebinger *et al.*, *Cell* 2019 — [doi:10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006) |
| IDEAS code | <https://github.com/maddalenastn/IDEAS> |
| IDEAS paper | Staiano *et al.*, *bioRxiv* 2025 — [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
| SCENT code (used here) | <https://github.com/Zaoqu-Liu/SCENT> |
| SCENT code (original) | <https://github.com/aet21/SCENT> |
| CCAT / signalling entropy paper | Teschendorff & Enver, *Nat Commun* 2017;8:15599 — [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599) |
| Waddington-OT | Schiebinger *et al.* 2019, above |

Every link above was resolved and checked, not written from memory.
