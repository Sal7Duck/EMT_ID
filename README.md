# Intrinsic dimension and signalling entropy in EMT

Code and data for the intrinsic dimension (ID) and signalling entropy analyses.

```
README.md                  this file — how every score was computed
environment.yml            plotting environment — the figure notebooks
environment-ideas.yml      IDEAS environment — the ID runs only

Figure_mouse_reprogramming/   13 panels — GSE115943
Figure_scc_12tumour/           5 panels — GSE288137
Figure_scc_multiome/           3 panels — GSE288220, RNA modality
Figure_hpcs/                   2 panels — GSE154989
Figure_transdiff_reprog/       2 panels — GSE112004

wot/                       Waddington-OT fate probabilities and trajectory weights
scores/                    the ID, entropy and ssGSEA score files
src/                       IDEAS, SCENT, the run scripts, the original analysis notebooks
```

Each `Figure_*/` folder is self-contained: a notebook, the input CSV it reads, a README
recording how every score was generated, and the rendered panels. **Open a notebook and run
it top to bottom** — no `.h5ad` is needed.

The pipeline runs in two stages, in two environments.

**Stage 1 — intrinsic dimension.** raw-count objects → `src/run_scripts/` → `scores/`.
Run in the `ideas` environment. Slow: the mouse case is 181,124 cells at k = 600.

**Stage 2 — figures.** `scores/` → `Figure_*/` → `Figure_*/figures/` and `source_data/`.
Run in the `emt-id` environment. Fast, and the only stage you need to make a figure.

The notebooks never call IDEAS. Each one shows the code that produced its scores, verbatim,
in section 1, and the `compute_ID` source in section 2 — so the calculation is documented
where you read the figure, without the notebook pretending it can re-run it.




---

# Methods

## Intrinsic dimension

ID was estimated with **IDEAS**, <https://github.com/maddalenastn/IDEAS>. IDEAS normalises
each cell to 10⁶ counts internally, so every input object holds **raw counts, not
log-transformed values**.

> Staiano M, Cirone N, Biondo M, Osella M, Scialdone A. Intrinsic dimensionality of
> single-cell transcriptomic data reveals potency landscapes during cell reprogramming.
> *bioRxiv* 2025.07.21.665922. [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922)

Two estimators are used.

**Per-cell (`local_2nn`)** — for each cell, the two-nearest-neighbour estimator of Facco et al.
is applied to that cell and its *k* nearest neighbours, giving one ID value per cell.

**Cluster-level (`2nn`)** — one ID per random subsample of a labelled group. Groups are
compared at equal *n*, drawn repeatedly, so that differences in group size cannot drive the
estimate.

`id_score=False` throughout, so IDEAS returns the raw ID. Any 0–1 values in the figures come
from min–max scaling applied afterwards in the notebook.

### Per dataset

| Dataset | Input object | Estimator | k / sampling | Output |
|---|---|---|---|---|
| Mouse reprogramming (GSE115943) | `IDEAS_input_mouse_reprogramming.h5ad` | `local_2nn` | k = 600 | `IDEAS_local_2NN_mouse_reprogramming.csv` |
| SCC, 12 tumours (GSE288137) | `emt_scc_ideas_input.h5ad` | `local_2nn` | k = 600 | `IDEAS_local_2NN_SCC_600.csv` |
| SCC, 12 tumours — Figure 3a | `emt_scc_ideas_input.h5ad` | `2nn` | 100 × 2,886 cells per state | `IDEAS_global_IDscore_100runs.csv` |
| SCC multiome (GSE288220) | `Multiome_scc_IDEAS_input.h5ad` | `local_2nn` | k = 600 | `IDEAS_local_2NN_multiome_SCC_600.csv` |
| HPCS (GSE154989) | `ideas_input_GSE154989_core_timecourse.h5ad` | `local_2nn` | k = 100 | `IDEAS_local_2NN_GSE154989_k100.csv` |

**Per-cell ID.**

```python
from computeID.functions_IDEAS import compute_ID

id_result = compute_ID(
    adata,                      # raw counts
    method="local_2nn",
    roots=None,
    n_neighbors=600,            # 100 for HPCS
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,
)
```

**Cluster-level ID (Figure 3a).** The four EMT states differ in size from 3,608 to 9,985
cells. Because ID estimates depend on sample size, states were compared at equal *n*:
2,886 cells — 80 % of the smallest state — drawn at random, 100 times per state.

```python
id_reps = compute_ID(
    adata,                      # raw counts
    group="celltype_annotation",
    method="2nn",
    sample_size=2886,
    n_samples=100,
    normalization=True,
    norm_sum=1e6,
    full_output=True,
    id_score=False,
)                               # -> 100 x 4 table
```

k = 600 for every dataset except HPCS, where k = 100. Local ID increases with k, so ID
values are compared within a dataset and not across datasets.

## Signalling entropy

Signalling entropy is the **CCAT** fast approximation, computed with `SCENT::CompCCAT` on
the log-normalised expression matrix against a protein–protein interaction network
(`net13Jun12`). The wrapper is in `src/scent/`. Provenance — version, commit, network, gene
count and orthology mapping — is written into `adata.uns["SCENT_CCAT"]` of each object.

> Method: Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
> differentiation potency from a cell's transcriptome. *Nat Commun* 2017;8:15599.
> [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)

Implementation used: <https://github.com/Zaoqu-Liu/SCENT>, an Rcpp-optimised fork of the
authors' original package <https://github.com/aet21/SCENT>.

## The score files

Two kinds of file live in `scores/`.

**The raw IDEAS output**, exactly as `compute_ID` returned it — this is what the run scripts
wrote, and it is what section 3 of each notebook loads:

| File | Rows |
|---|---|
| `mouse_reprogramming/IDEAS_local_2NN_mouse_reprogramming.csv` | 181,124 barcodes |
| `scc_12tumour/IDEAS_local_2NN_SCC_600.csv` | 23,959 barcodes |
| `scc_12tumour/IDEAS_global_IDscore_100runs.csv` | 100 subsamples x 4 states |
| `scc_multiome/IDEAS_local_2NN_multiome_SCC_600.csv` | 14,302 barcodes |
| `hpcs/IDEAS_local_2NN_GSE154989_k100.csv` | 3,891 barcodes |

**The cell tables**, `scores/<dataset>/cell_table.csv.gz` — one row per barcode, holding the
ID score together with the signalling entropy, the other scores, and the annotations
(cell type, EMT state, tumour, mouse, time point, WOT weights). Built once by
`src/export_cell_tables.py` from the annotated `.h5ad` objects.

| Table | Cells | Key columns |
|---|---|---|
| `mouse_reprogramming/cell_table.csv.gz` | 181,124 | `timepoint`, `IDEAS_score`, `CCAT_score`, `KS_score`, `GS76_score`, `score_ks_epithelial`, `score_ks_mesenchymal`, `score_pluripotency`, `w_*`, `P_*` |
| `scc_12tumour/cell_table.csv.gz` | 23,959 | `cluster_annotation`, `cluster`, `sample`, `assay`, `ID_score`, `CCAT_score`, `transcriptome_entropy`, `HPCS_ssGSEA` |
| `scc_multiome/cell_table.csv.gz` | 14,302 | `emt_state`, `ID_score` |
| `hpcs/cell_table.csv.gz` | 3,891 | `cell_type`, `lineage`, `mouseID`, `genotype`, `week`, `IDEAS_local_ID_100` |

The index is the barcode exactly as it appears in the `.h5ad`, so a table joins back onto the
object with a plain index merge. To rebuild them:

```bash
python src/export_cell_tables.py --data-dir /path/to/h5ads
```

## Figure notebooks

One notebook per figure page. Each one shows the verbatim code that produced its ID scores,
the `compute_ID` source, the score file, the cell table, then the panels and the statistics.

**The notebooks read only the CSVs.** No `.h5ad` is needed to make a figure — the 8.6 GB
mouse object is used once, by `export_cell_tables.py`, and never again.

| Notebook | Reproduces |
|---|---|
| `Figure_mouse_reprogramming/` | ID across time, WOT trajectories, six correlation panels, entropy |
| `Figure_scc_12tumour/` | cluster-level ID across 4 EMT states, KS scatter, HPCS programme |
| `Figure_scc_multiome/` | cluster-level ID across 5 EMT states, KS scatter |
| `Figure_hpcs/` | ID across 12 lung tumour cell types, HPCS effect size |
| `Figure_transdiff_reprog/` | ID and entropy along transdifferentiation vs reprogramming |

The right-hand KS-Epi vs KS-Mes scatters on pages 3a and 3b need the ssGSEA KS scores, which
were computed inside the analysis notebooks and never saved. `src/compute_ks_ssgsea.py`
recomputes them from the expression matrix:

```bash
python src/compute_ks_ssgsea.py <object.h5ad> scores/<dataset>/KS_ssGSEA.csv
```

Done for the 12-tumour SCC object — `scores/scc_12tumour/KS_ssGSEA.csv`, which reproduces the
published r = -0.627. The multiome object is not on this machine yet, so page 3b's right panel
is still missing.

## Statistics

**Group comparisons.** Kruskal–Wallis across groups, then Dunn's post-hoc with pooled ranks,
tie correction and Benjamini–Hochberg FDR. Effect sizes are reported as median differences
and rank-biserial correlation. For the reprogramming time course, where the comparison is
across ordered time bins, one-way ANOVA with Tukey HSD is used instead.

**Correlations.** Pearson *r* on cells with both scores present, with the ordinary
least-squares line drawn through the same cells. At these sample sizes the p-value underflows
the double-precision floor, and is reported as `p < 1e-300`.

**Waddington-OT weighting.** Trajectory means are weighted by the WOT ancestor probability
`w`, renormalised within each day over cells with non-zero weight. The effective sample size
is `n_eff = (Σw)² / Σw²`.

---

## References

| | |
|---|---|
| Mouse reprogramming | [GSE115943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115943) · Schiebinger *et al.*, *Cell* 2019;176:928–943.e22 · [doi:10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006) |
| SCC, 12 tumours | [GSE288137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288137) |
| SCC multiome | [GSE288220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288220) |
| HPCS lung | [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989) |
| IDEAS | <https://github.com/maddalenastn/IDEAS> · Staiano *et al.*, *bioRxiv* 2025 · [doi:10.1101/2025.07.21.665922](https://doi.org/10.1101/2025.07.21.665922) |
| SCENT / CCAT | <https://github.com/Zaoqu-Liu/SCENT> (used) · <https://github.com/aet21/SCENT> (original) · Teschendorff & Enver, *Nat Commun* 2017;8:15599 · [doi:10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599) |
