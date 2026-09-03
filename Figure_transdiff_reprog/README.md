# Figure_transdiff_reprog

Intrinsic dimension and signalling entropy along two trajectories out of the same starting
cell type — one transdifferentiating, one reprogramming.

```
Figure_transdiff_reprog/
  README.md
  Figure_transdiff_reprog.ipynb   the notebook — two panels
  Input_transdiff_reprog.csv      2,688 barcodes
  figures/                        600 dpi PNG + vector PDF
  source_data/                    the numbers behind each panel
```

Environment: `../environment.yml` (`emt-id`).

---

## Where the data came from

**[GSE112004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112004)** — "Single cell
RNA-seq identifies the origins of heterogeneity in efficient cell transdifferentiation and
reprogramming". Mouse pre-B cells from bone marrow, 3,840 cells, Illumina HiSeq 2500.

> Francesconi M, Di Stefano B, Berenguer C, de Andrés-Aguayo L, Plana-Carmona M,
> Mendez-Lago M, Guillaumet-Adkins A, Rodriguez-Esteban G, Gut M, Gut IG, Heyn H, Lehner B,
> Graf T. *eLife* 2019;8:e41627.
> [doi:10.7554/eLife.41627](https://doi.org/10.7554/eLife.41627) · PMID 30860479

Source scores: `../scores/gse112004/GSE112004_scores.csv`.

### Why this dataset is the right control

The same pre-B cell is pushed down two routes:

| Trajectory | Induction | Treatment | Ends at | Time points |
|---|---|---|---|---|
| **iMac** | C/EBPα | transdifferentiation | induced macrophage | 42h, 66h, 114h |
| **iPSC** | OSKM | reprogramming | induced pluripotent cell | day 2, 4, 6, 8 |

One route ends differentiated, the other pluripotent, from a common origin. If the intrinsic
dimension tracks potency it should fall along the first and rise along the second — which is
what the figure shows.

### What is dropped, and why

The three earliest time points — **0h, 6h, 18h**, 1,152 cells — precede the fate split and
carry no trajectory label in the source file. They are excluded; **2,688 cells** remain,
384 at each of the 7 time points.

---

## How each score was generated

### `ID_score` — intrinsic dimension

**Per-cell** IDEAS, `method="local_2nn"`, **k = 100**. One ID per cell, estimated from that
cell's own 100-nearest-neighbour neighbourhood.

**IDEAS** — <https://github.com/maddalenastn/IDEAS>.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922.
> [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)

```python
id_local = compute_ID(
    adata,                      # 3840 x 28874, raw UMI counts
    group="all",
    method="local_2nn",
    id_score=True,
    n_neighbors=100,
    normalization=True,
    norm_sum=int(1e6),
)
```

**`id_score=True` here**, unlike the other figures in this repository, which use `False` and
scale afterwards. IDEAS therefore returned values already scaled 0–1 across all 3,840 cells.
They are rescaled over the 2,688 plotted cells so each panel spans the full range.
Both steps are linear, so no statistic changes.

**A global 2NN ID is also in the source file** (`ID_global_2nn`: one value per time point,
repeated across that group's cells, 30 subsamples). It is **not used in this figure** — the
panels are the per-cell local ID.

### `signalling_entropy` — CCAT

`SCENT::CompCCAT` on the log-normalised matrix against a protein–protein interaction network.

> Teschendorff AE, Enver T. *Nat Commun* 2017;8:15599.
> [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)

Code used: <https://github.com/Zaoqu-Liu/SCENT>, an Rcpp-optimised fork of
<https://github.com/aet21/SCENT>. Wrapper in `../src/scent/`.

---

## Scaling

Both scores are min–max scaled 0–1 over the 2,688 plotted cells. Raw values are kept in
`ID_score_raw` and `signalling_entropy_raw`. Min–max is linear, so Kruskal-Wallis *H*,
Dunn *q* and Mann-Whitney *U* are unchanged.

---

## Statistics

**Within each trajectory** — Kruskal-Wallis across that trajectory's time points, then Dunn
with pooled ranks, tie correction and BH-FDR.

**Between trajectories** — Mann-Whitney on the endpoints (iMac 114h vs iPSC day 8), with
rank-biserial correlation as the effect size.

### What the tests say

| | iMac | iPSC |
|---|---|---|
| ID, Kruskal-Wallis | H = 433.6, p = 7e-95 | H = 634.2, p = 4e-137 |
| Signalling entropy, Kruskal-Wallis | H = 81.1, p = 2e-18 | H = 702.2, p = 7e-152 |

**The two routes separate completely at the endpoint.** ID: median 0.180 (iMac 114h) vs 0.481
(iPSC day 8), rank-biserial −0.974 — near-total separation of the two distributions.
Signalling entropy: 0.417 vs 0.714, rank-biserial −0.824.

Both measures move the same way: **down** along transdifferentiation, **up** along
reprogramming. That the two independent measures agree, on the same cells from the same
origin, is the point of the panel.

One nuance worth stating: the ID rises from 42h to 66h before falling at 114h, so the iMac
decline is not monotonic — 42h vs 66h is significant in the opposite direction.

---

## References

| | |
|---|---|
| Dataset | [GSE112004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112004) |
| Dataset paper | Francesconi *et al.*, *eLife* 2019;8:e41627 — [doi:10.7554/eLife.41627](https://doi.org/10.7554/eLife.41627) |
| IDEAS | <https://github.com/maddalenastn/IDEAS> · [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
| SCENT / CCAT | <https://github.com/Zaoqu-Liu/SCENT> · <https://github.com/aet21/SCENT> · [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599) |
