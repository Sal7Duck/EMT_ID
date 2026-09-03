# Intrinsic dimension runs

These scripts compute the ID scores. They are run **separately from the figure notebooks**, in
the IDEAS environment:

```bash
conda env create -f ../../environment-ideas.yml
conda activate ideas
python run_mouse_ideas.py
```

Each script reads an input object from `data/inputs/`, calls `compute_ID`, and writes a score
file into `scores/`. The figure notebooks read those score files and never call IDEAS.

| Script | Input | Method | Output |
|---|---|---|---|
| `run_mouse_ideas.py` | `IDEAS_input_mouse_reprogramming.h5ad` | `local_2nn`, k = 600 | `scores/mouse_reprogramming/IDEAS_local_2NN_mouse_reprogramming.csv` |
| `run_scc_ideas.py` | `emt_scc_ideas_input.h5ad` | `local_2nn`, k = 600 | `scores/scc_12tumour/IDEAS_local_2NN_SCC_600.csv` |
| `run_scc_cluster_ideas.py` | `emt_scc_ideas_input.h5ad` | `2nn`, 100 × 2,886 cells | `scores/scc_12tumour/IDEAS_global_IDscore_100runs.csv` |
| `run_ideas_multiome.py` | `Multiome_scc_IDEAS_input.h5ad` | `local_2nn`, k = 600 | `scores/scc_multiome/IDEAS_local_2NN_multiome_SCC_600.csv` |
| `run_hpcs_ideas.py` | `ideas_input_GSE154989_core_timecourse.h5ad` | `local_2nn`, k = 100 | `scores/hpcs/IDEAS_local_2NN_GSE154989_k100.csv` |

`run_ideas.py` and `run_id.py` are earlier runs kept for the record; they did not produce the
deposited scores.

## The two methods

**`local_2nn`** returns one ID per cell, so the result merges into `adata.obs` by cell
identifier. Used for figures 1, 3b and 3c.

**`2nn`** returns one ID per random subsample of a group — a `n_samples × n_groups` table with
no cell identity. It cannot be merged into `obs`, and is used as its own table. Figure 3a.

Paths inside the scripts are the absolute cluster paths they were run with. Edit them, or run
from a directory laid out the same way.
