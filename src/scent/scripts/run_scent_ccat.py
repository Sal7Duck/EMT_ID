import os
import subprocess
import tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import mygene


PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

R_BIN = "/home/csb/miniconda3/envs/scent-r/bin/Rscript"

R_SCRIPT = os.path.join(
    PROJECT_ROOT,
    "scripts",
    "run_ccat.R"
)

PPI_FILE = os.path.join(
    PROJECT_ROOT,
    "SCENT",
    "data",
    "net13Jun12.rda"
)

ORTHOLOGY_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "reference",
    "mouse_human_1to1.csv"
)


def build_mouse_human_mapping(mouse_genes):
    """
    Build a strict 1:1 mouse-symbol -> human-Entrez mapping
    using MyGene/HomoloGene.

    The resulting mapping is saved locally so that future datasets
    can reuse exactly the same reference mapping.
    """

    mouse_genes = pd.Series(
        mouse_genes, dtype="string"
    ).dropna().unique().tolist()

    mg = mygene.MyGeneInfo()

    results = mg.querymany(
        mouse_genes,
        scopes="symbol",
        fields=["symbol", "entrezgene", "homologene"],
        species="mouse",
        as_dataframe=False,
        returnall=False,
        verbose=False
    )

    candidates = []

    for r in results:

        mouse_symbol = r.get("symbol")

        if mouse_symbol is None:
            mouse_symbol = r.get("query")

        homologene = r.get("homologene")

        if not homologene:
            continue

        for pair in homologene.get("genes", []):

            if len(pair) != 2:
                continue

            taxid, gene_id = pair

            # Human taxon
            if str(taxid) == "9606":
                candidates.append(
                    {
                        "mouse_symbol": str(mouse_symbol),
                        "human_entrez": str(gene_id)
                    }
                )

    candidates = pd.DataFrame(candidates)

    if candidates.empty:
        raise RuntimeError(
            "No mouse-human orthologs were found."
        )

    # Remove exact duplicates
    candidates = candidates.drop_duplicates()

    # Keep only mouse genes with exactly ONE human candidate
    mouse_counts = (
        candidates.groupby("mouse_symbol")
        ["human_entrez"]
        .nunique()
    )

    candidates = candidates[
        candidates["mouse_symbol"].isin(
            mouse_counts[mouse_counts == 1].index
        )
    ]

    # Keep only human Entrez IDs with exactly ONE mouse gene
    human_counts = (
        candidates.groupby("human_entrez")
        ["mouse_symbol"]
        .nunique()
    )

    mapping = candidates[
        candidates["human_entrez"].isin(
            human_counts[human_counts == 1].index
        )
    ].copy()

    mapping = mapping.sort_values(
        ["mouse_symbol", "human_entrez"]
    ).reset_index(drop=True)

    os.makedirs(
        os.path.dirname(ORTHOLOGY_FILE),
        exist_ok=True
    )

    mapping.to_csv(
        ORTHOLOGY_FILE,
        index=False
    )

    print(
        f"Saved mouse-human mapping: {ORTHOLOGY_FILE}"
    )

    print(
        f"Strict 1:1 mappings: {len(mapping)}"
    )

    return mapping


def load_or_build_mapping(mouse_genes):
    """
    Reuse the saved mapping if it exists.
    Otherwise construct it.
    """

    if os.path.exists(ORTHOLOGY_FILE):

        mapping = pd.read_csv(
            ORTHOLOGY_FILE,
            dtype=str
        )

        print(
            f"Using existing mapping: {ORTHOLOGY_FILE}"
        )

    else:

        mapping = build_mouse_human_mapping(
            mouse_genes
        )

    return mapping


def run_scent_ccat(
    adata,
    species="mouse",
    counts_layer="counts",
    chunk_size=1000
):
    """
    Calculate SCENT CCAT scores using the official
    SCENT CompCCAT implementation.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.

    species : str
        "mouse" or "human".

    counts_layer : str
        Layer containing raw counts.

    chunk_size : int
        Number of cells sent to R per iteration.

    Returns
    -------
    AnnData
        Copy of adata with CCAT_score in obs.
    """

    if counts_layer not in adata.layers:
        raise ValueError(
            f"Layer '{counts_layer}' not found."
        )

    species = species.lower()

    if species not in ["mouse", "human"]:
        raise ValueError(
            "species must be 'mouse' or 'human'"
        )

    # --------------------------------------------------
    # Gene mapping
    # --------------------------------------------------

    if species == "mouse":

        mapping = load_or_build_mapping(
            adata.var_names
        )

        mapping = mapping[
            mapping["mouse_symbol"].isin(
                adata.var_names
            )
        ].copy()

        print(
            f"AnnData genes: {adata.n_vars}"
        )

        print(
            f"Genes with strict 1:1 mapping: "
            f"{len(mapping)}"
        )

        if len(mapping) < 5000:
            raise RuntimeError(
                "Fewer than 5000 mouse genes have "
                "strict 1:1 human orthologs."
            )

        # Map AnnData genes to human Entrez IDs
        gene_to_entrez = dict(
            zip(
                mapping["mouse_symbol"],
                mapping["human_entrez"]
            )
        )

        selected_mouse = [
            g for g in adata.var_names
            if g in gene_to_entrez
        ]

        human_entrez = [
            gene_to_entrez[g]
            for g in selected_mouse
        ]

    else:

        selected_mouse = list(
            adata.var_names.astype(str)
        )

        human_entrez = selected_mouse

    # --------------------------------------------------
    # Check duplicate human IDs
    # --------------------------------------------------

    tmp = pd.DataFrame(
        {
            "gene": selected_mouse,
            "entrez": human_entrez
        }
    )

    duplicated = tmp["entrez"].duplicated(
        keep=False
    )

    if duplicated.any():

        # Keep first occurrence of each human Entrez ID
        keep = ~tmp["entrez"].duplicated(
            keep="first"
        )

        tmp = tmp.loc[keep].copy()

    selected_mouse = tmp["gene"].tolist()
    human_entrez = tmp["entrez"].tolist()

    print(
        f"Genes passed toward SCENT: "
        f"{len(human_entrez)}"
    )

    # --------------------------------------------------
    # Prepare output
    # --------------------------------------------------

    scores = pd.Series(
        np.nan,
        index=adata.obs_names,
        dtype=float
    )

    # --------------------------------------------------
    # Process cells in chunks
    # --------------------------------------------------

    for start in range(
        0,
        adata.n_obs,
        chunk_size
    ):

        end = min(
            start + chunk_size,
            adata.n_obs
        )

        print(
            f"CCAT: cells "
            f"{start:,} - {end:,} "
            f"of {adata.n_obs:,}"
        )

        X = adata[
            start:end,
            selected_mouse
        ].layers[counts_layer]

        if hasattr(X, "toarray"):
            X = X.toarray()

        X = np.asarray(
            X,
            dtype=np.float64
        )

        # SCENT expects genes x cells
        expr = pd.DataFrame(
            X.T,
            index=human_entrez,
            columns=adata.obs_names[
                start:end
            ]
        )

        with tempfile.TemporaryDirectory() as tmpdir:

            expression_file = os.path.join(
                tmpdir,
                "expression.tsv"
            )

            output_file = os.path.join(
                tmpdir,
                "ccat.tsv"
            )

            expr.to_csv(
                expression_file,
                sep="\t"
            )

            command = [
                R_BIN,
                R_SCRIPT,
                expression_file,
                PPI_FILE,
                output_file
            ]

            subprocess.run(
                command,
                check=True
            )

            result = pd.read_csv(
                output_file,
                sep="\t"
            )

        result = result.set_index(
            "cell_id"
        )

        scores.loc[
            result.index
        ] = result["CCAT_score"]

    # --------------------------------------------------
    # Add to AnnData
    # --------------------------------------------------

    adata.obs["CCAT_score"] = scores

    # --------------------------------------------------
    # Record provenance
    # --------------------------------------------------

    adata.uns["SCENT_CCAT"] = {
        "SCENT_version": "2.0.0",
        "SCENT_commit":
            "8a3a80953befa50b4320249d28e1ce0afb6c5e06",
        "method": "CompCCAT",
        "ppi": "net13Jun12",
        "ppi_file": PPI_FILE,
        "species": species,
        "counts_layer": counts_layer,
        "orthology": (
            "Mouse symbol -> human Entrez "
            "via MyGene/HomoloGene"
            if species == "mouse"
            else "Not applicable"
        ),
        "mapping_type": (
            "strict 1:1"
            if species == "mouse"
            else "human Entrez"
        ),
        "n_genes_used": len(human_entrez),
        "chunk_size": chunk_size
    }

    return adata
