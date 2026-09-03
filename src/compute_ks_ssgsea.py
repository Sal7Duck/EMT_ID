#!/usr/bin/env python3
"""
KS epithelial / KS mesenchymal ssGSEA scores for the SCC datasets.

These two scores are what the right-hand panels of figures 3a and 3b plot. They were
computed inside the analysis notebooks and never written to the deposited objects, so this
script recomputes them from the expression matrix and saves them as CSVs keyed by barcode.

Verbatim from `id_entropy_12_emt_scc.ipynb` cell 80 and `Multiome_cluster_id.ipynb` cell 14:
the same two mouse signatures, the same gseapy call, the same chunking.

`sample_norm_method="rank"` ranks the genes within each cell, so it does not matter whether
X holds raw counts or log-normalised values — log1p and total-count scaling are monotone
within a cell and leave the ranks unchanged.

    python src/compute_ks_ssgsea.py <input.h5ad> <output.csv>

Needs gseapy:  pip install gseapy
"""

import os
import sys

import anndata as ad
import gseapy as gp
import numpy as np
import pandas as pd
from scipy.sparse import issparse

CHUNK_SIZE = 500      # cells per ssGSEA call; keeps peak memory near 1 GB
THREADS = 4

# Mouse KS epithelial gene signature: 164 genes
KS_EPI_MOUSE = """
Cdh1 Agr2 Epcam Krt19 Rab25 Tacstd2 Galnt3 Fxyd3 Spint2 Tmem30b Scnn1a St14
Esrp1 S100a14 Cldn7 Erbb3 Rbm47 Spint1 Elf3 Cldn13 Prss8 Sh3yl1 Ehf Lcn2 Jup
Vamp8 Krt8 Inava Krt7 Dsp Cds1 Itgb4 Tmprss4 Lsr Sorl1 Grhl2 Ppl AA986860
Tspan1 Map7 Slpi Tox3 Arhgap8 F11r Lad1 Gpx2 Ctsh Adgrg1 Fa2h Klf5 Areg Krt18
Scel Cdh3 Ugt1a1 Mpzl2 Crybg1 Ovol2 Llgl2 Esrp2 Myo5c Ddr1 Vgll1 Irf6 Sfn
Tspan13 Kcnk1 Myo1d Pkp3 Itgb6 Ly75 Mapk13 Ttc39a Elmo3 Ceacam1 Dtx4 Erbb2
Rab11fip1 Atp2c2 Mst1r Ap1m2 Tgfa Myo6 Ptk6 Oas1g Fbp1 Aqp3 Cblc Epha1 Bspry
Eps8l1 Grb7 0610040J01Rik Klk6 Tjp2 Pls1 Dennd2d Eps8l2 Il20ra Hes1 Il1rn
Exph5 Arhgdib Misp Camk2n1 Hpgd Syngr2 Perp Mansc1 Dsc2 Pof1b Serinc5 Bik
Anxa9 Mall Epn3 Stap2 Foxa1 Pycard Slc37a1 Ank3 Tspan15 Hnmt Abcc3 Sdc1 Ckmt1
Tob1 B3gnt3 Tmc6 Cd9 Adap1 Atp1b1 Shank2 Cyb561 Ermp1 Rab20 Myh14 Capn1
Aldh3b2 Trim31 Arap2 Ssh3 Ica1 Arhgef5 Alox5 Rhod Tmprss2 Mtus1 Cyp4f18
Ppfibp2 Rabgap1l Plxnb2 Mgst2 Evpl Cd46 Krt15 Cnksr1 Blnk Comt Anxa4 Tnfsf13
Ocln Nherf1 Xbp1
""".split()

# Mouse KS mesenchymal gene signature: 48 genes
KS_MES_MOUSE = """
Vim Zeb1 Emp3 Sacs Axl Loxl2 Sparc Fhl4 Fermt2 Tuba1a Tmem158 Cald1 Lgals1
Pmp22 Msn Glyr1 Map1b Ap1s2 Gja1 Dennd5a Fam216a Tpm2 Tubb6 Srpx Ank2 Chn1
Sh2b3 P3h1 Etv1 Sobp Akap12 Tgfb1i1 Serpine1 Soat1 Col5a2 Lhfpl6 Cep170 Popdc3
Trpc1 Poglut2 Myl9 Bag2 Fstl1 Mxra7 Gfpt2 Reck Tmeff1 Cavin1
""".split()


def main(input_h5ad, output_csv):
    # backed="r": the expression matrix stays on disk and is read one chunk at a
    # time. A 1.7 GB object does not fit in memory alongside the dense chunk.
    adata = ad.read_h5ad(input_h5ad, backed="r")
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes")

    available = set(adata.var_names.astype(str))
    epi_found = [g for g in KS_EPI_MOUSE if g in available]
    mes_found = [g for g in KS_MES_MOUSE if g in available]
    print(f"KS-Epi genes found: {len(epi_found)}/{len(KS_EPI_MOUSE)}")
    print(f"KS-Mes genes found: {len(mes_found)}/{len(KS_MES_MOUSE)}")

    if len(epi_found) < 20 or len(mes_found) < 5:
        raise ValueError(
            "Too few signature genes matched. var_names must be mouse gene symbols "
            "(Cdh1, Epcam, Vim), not Ensembl IDs.")

    gene_sets = {"KS_Epi": KS_EPI_MOUSE, "KS_Mes": KS_MES_MOUSE}

    # Checkpoint each chunk next to the output, so an interrupted run resumes
    # instead of starting over. 24 chunks at ~30 s each.
    ckpt_dir = output_csv + ".chunks"
    os.makedirs(ckpt_dir, exist_ok=True)
    all_scores = []

    for start in range(0, adata.n_obs, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, adata.n_obs)

        ckpt = os.path.join(ckpt_dir, f"{start:07d}.csv")
        if os.path.exists(ckpt):
            cached = pd.read_csv(ckpt, index_col=0)
            # A checkpoint written under a different CHUNK_SIZE has the wrong length;
            # recompute rather than trust it.
            if len(cached) == stop - start:
                all_scores.append(cached)
                print(f"cells {start + 1:,}-{stop:,}: from checkpoint", flush=True)
                continue
            print(f"cells {start + 1:,}-{stop:,}: stale checkpoint "
                  f"({len(cached)} rows), recomputing", flush=True)
        X = adata.X[start:stop]
        if hasattr(X, "to_memory"):          # backed sparse view
            X = X.to_memory()
        X = X.T.toarray() if issparse(X) else np.asarray(X).T

        expression = pd.DataFrame(X,
                                  index=adata.var_names.astype(str),
                                  columns=adata.obs_names[start:stop].astype(str))
        expression = expression.groupby(level=0, sort=False).mean()

        result = gp.ssgsea(data=expression, gene_sets=gene_sets,
                           sample_norm_method="rank", permutation_num=0,
                           min_size=5, max_size=500, outdir=None,
                           no_plot=True, threads=THREADS)

        chunk = (result.res2d.pivot(index="Name", columns="Term", values="ES")
                 .reindex(adata.obs_names[start:stop].astype(str)))
        chunk.to_csv(ckpt)
        all_scores.append(chunk)
        print(f"scored cells {start + 1:,}-{stop:,} of {adata.n_obs:,}", flush=True)

    # Concatenate in cell order. No final reindex: barcodes repeat across tumours in
    # the SCC object, and reindexing on a duplicated axis raises.
    scores = pd.concat(all_scores)
    assert len(scores) == adata.n_obs, f"{len(scores)} rows for {adata.n_obs} cells"

    out = pd.DataFrame({"KS_Epi_ssGSEA": scores["KS_Epi"].astype(float).to_numpy(),
                        "KS_Mes_ssGSEA": scores["KS_Mes"].astype(float).to_numpy()},
                       index=adata.obs_names.astype(str))
    out.index.name = "barcode"
    out.to_csv(output_csv)
    print(f"\nwritten: {output_csv}")
    print(out.describe().round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
