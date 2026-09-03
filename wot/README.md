# Waddington-OT

How the fate probabilities and trajectory weights used in
`../Figure_mouse_reprogramming/` were computed.

```
wot/
  README.md                    this file
  fate_probabilities.csv.gz    P — 181,124 cells x 5 fates
  trajectory_weights.csv.gz    w — 181,124 cells x 4 fates
```

Both are also carried in `../Figure_mouse_reprogramming/Input_mouse_reprogramming.csv` as
the `P_*` and `w_*` columns, so the figure needs nothing from this folder. They are here as
the standalone record.

---

## Method

Trajectories were computed with **Waddington-OT** (`wot`), which fits entropically
regularised, unbalanced optimal transport maps between consecutive time points and composes
them to follow populations through the time course.

> Schiebinger G, Shu J, Tabaka M, Cleary B, Subramanian V, Solomon A, Gould J, Liu S, Lin S,
> Berube P, Lee L, Chen J, Brumbaugh J, Rigollet P, Hochedlinger K, Jaenisch R, Regev A,
> Lander ES. Optimal-transport analysis of single-cell gene expression identifies
> developmental trajectories in reprogramming. *Cell* 2019;176:928–943.e22.
> [doi:10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006)
>
> Software: <https://broadinstitute.github.io/wot/>

Maps were fitted with `wot.ot.OTModel` on the 3,000 highly variable genes, over all 181,124
cells and 39 time points from day 0 to day 18, giving **38 transport maps** between adjacent
time points.

```python
ot_model = wot.ot.OTModel(
    adata[:, adata.var["highly_variable"]],
    day_field="day",
    epsilon=0.05,
    lambda1=1,
    lambda2=50,
    growth_iters=3,
    local_pca=30,
)
ot_model.compute_all_transport_maps(tmap_out="tmaps/serum")

tmap_model  = wot.tmap.TransportMapModel.from_directory("tmaps/serum")
populations = tmap_model.population_from_cell_sets(cell_sets, at_time=18.0)

P = tmap_model.fates(populations)          # fate_probabilities.csv.gz
w = tmap_model.trajectories(populations)   # trajectory_weights.csv.gz
```

### Parameters

| Parameter | Value | What it controls |
|---|---|---|
| `epsilon` | 0.05 | entropic regularisation — how diffuse each map is |
| `lambda1` | 1 | unbalancedness on the source distribution |
| `lambda2` | 50 | unbalancedness on the target distribution |
| `growth_iters` | 3 | rounds of growth-rate refinement |
| `local_pca` | 30 | dimensions of the local cost metric |
| highly variable genes | 3,000 | the space the maps are fitted in |
| time points | 39 | day 0 to day 18 |
| transport maps | 38 | one per adjacent pair, 7.08 GB total |

Growth rates use the standard WOT birth–death form driven by the Proliferation and Apoptosis
signature scores.

### Gene sets

**All gene sets are taken from the original paper** — Schiebinger *et al.* 2019, Table S2
(`mmc2-2.xlsx`), unmodified. Six sets are used: Proliferation and Apoptosis for the growth
rates, and four identity signatures for the terminal fates. They are scored with
`wot.score_gene_sets(method="mean_z_score", max_z_score=5)`.

---

## The four fates

Terminal populations are the **top 500 cells at day 18** by each identity signature.

| Fate | Gene set (Table S2) | Cells | Score cutoff | Overlap with other sets |
|---|---|---|---|---|
| iPSC | Pluripotency | 500 | 0.837 | 3 (Neural) |
| Stromal | MEF identity | 500 | 1.353 | 28 (Neural) |
| Trophoblast | Trophoblast identity | 500 | 0.691 | none |
| Neural | Neural identity | 500 | 0.088 | 31 (iPSC 3, Stromal 28) |

2,000 assignments over **1,969 distinct cells** — 31 are claimed by two fates, all at the
Stromal/Neural boundary. The four sets are 15 % of the 13,386 cells at day 18.

This is the most consequential choice in the pipeline: every probability in the two files
below follows from these 2,000 cells.

---

## The two files

Both come from the same transport maps, normalised in opposite directions. They are **not**
interchangeable.

| File | Columns | Normalised | Answers |
|---|---|---|---|
| `fate_probabilities.csv.gz` — **P** | iPSC, Stromal, Trophoblast, Neural, Other | across fates, **per cell** → each row sums to 1 | *where is this cell going?* |
| `trajectory_weights.csv.gz` — **w** | iPSC, Stromal, Trophoblast, Neural | across cells, **per fate per day** → each column sums to 1 within a day | *which cells is this fate descended from?* |

`P` has the extra `Other` column, absorbing the mass that reaches none of the four
populations. Both files are indexed by cell barcode, matching the annotated object.
