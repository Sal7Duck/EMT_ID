# Original analysis notebooks

The notebooks the figure code was taken from, exactly as they were run — kept here so the
provenance quoted in each `Figure_*/` notebook points at a file you can actually open.

| Notebook | What it produced |
|---|---|
| `Figure_paper.ipynb` | the mouse reprogramming panels |
| `local_id.ipynb` | the 12-tumour SCC cluster-level ID (`method="2nn"`) |
| `id_entropy_12_emt_scc.ipynb` | the 12-tumour SCC KS ssGSEA scores and scatter |
| `Multiome_cluster_id.ipynb` | the multiome cluster-level ID and KS ssGSEA |

**These are the record, not the pipeline.** They carry absolute paths from the analysis
machine and will not run as-is. The runnable versions are `../run_scripts/` and
`../run_cluster_id.py`; the figures are rebuilt from `../../Figure_*/`.

Cell outputs other than printed text have been stripped — the embedded figures made these
files 5–8 MB each. The code is untouched.
