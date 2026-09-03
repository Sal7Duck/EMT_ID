# SCENT: Single Cell ENTropy for Differentiation Potency Estimation

<!-- badges: start -->
[![R-Universe](https://zaoqu-liu.r-universe.dev/badges/SCENT)](https://zaoqu-liu.r-universe.dev/SCENT)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/Zaoqu-Liu/SCENT)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub](https://img.shields.io/badge/GitHub-Zaoqu--Liu%2FSCENT-black?logo=github)](https://github.com/Zaoqu-Liu/SCENT)
<!-- badges: end -->

## Overview

**SCENT** is an R package for estimating differentiation potency of single cells from scRNA-Seq data using signaling entropy on protein-protein interaction (PPI) networks. The package implements the computational framework described in [Teschendorff & Enver (2017)](https://www.nature.com/articles/ncomms15599), with performance optimizations through C++ (Rcpp/RcppArmadillo).

📖 **Documentation**: https://zaoqu-liu.github.io/SCENT/

## Theoretical Background

### Signaling Entropy Rate (SR)

The differentiation potency of a cell can be quantified by the **signaling entropy rate**, which measures the "randomness" of information flow through the protein interaction network. The mathematical formulation is:

$$SR = \frac{1}{SR_{max}} \sum_{j} \pi_j S_j$$

where:

- **$\pi_j$**: Stationary distribution at gene $j$, computed as $\pi_j = \frac{x_j \cdot (\mathbf{A}\mathbf{x})_j}{Z}$
- **$S_j$**: Local entropy at gene $j$, defined as $S_j = -\sum_k P_{jk} \log P_{jk}$
- **$P_{jk}$**: Transition probability from gene $j$ to $k$, given by $P_{jk} = \frac{A_{jk} \cdot x_k}{\sum_l A_{jl} \cdot x_l}$
- **$SR_{max}$**: Maximum entropy rate, equal to $\log(\lambda_1)$ where $\lambda_1$ is the largest eigenvalue of $\mathbf{A}$

### Biological Interpretation

| Cell State | SR Value | Network Signaling |
|------------|----------|-------------------|
| Pluripotent (stem-like) | High | Uniform flow through multiple pathways |
| Differentiated | Low | Concentrated flow through specific pathways |

## Installation

### From R-Universe (Recommended)

```r
install.packages("SCENT", repos = c(
  "https://zaoqu-liu.r-universe.dev",
  "https://cloud.r-project.org"
))
```

### From GitHub

```r
# Install dependencies
install.packages(c("Rcpp", "RcppArmadillo", "Matrix", "igraph"))

# Install SCENT
devtools::install_github("Zaoqu-Liu/SCENT")
```

### System Requirements

| Platform | Requirement |
|----------|-------------|
| macOS | `xcode-select --install` |
| Windows | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) |
| Linux | Standard build tools (gcc, g++) |

## Quick Start

```r
library(SCENT)

# Load built-in PPI network (8,434 genes, 303,600 interactions)
data(net13Jun12.m)

# Load your expression matrix (genes × cells)
# For Seurat objects: exp <- GetAssayData(seurat_obj, layer = "data")

# Method 1: CCAT (fast approximation)
ccat_scores <- CompCCAT(exp, net13Jun12.m)

# Method 2: SR (full computation)
integ <- DoIntegPPI(exp, net13Jun12.m)
sr_result <- CompSRana(integ, local = TRUE)

# Access results
sr_scores <- sr_result$SR           # Global entropy rate per cell
local_entropy <- sr_result$locS     # Gene-level local entropy
stationary <- sr_result$inv         # Stationary distribution
```

## Methods Comparison

| Method | Algorithm | Speed | Accuracy | Use Case |
|--------|-----------|-------|----------|----------|
| **SR** | Full entropy computation | ~5s / 100 cells | Gold standard | Publication, validation |
| **CCAT** | Correlation with degree | ~0.1s / 100 cells | Approximation (r ≈ 0.78) | Screening, large datasets |

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [Quick Start](https://zaoqu-liu.github.io/SCENT/articles/quick-start.html) | Basic usage and workflow |
| [Algorithm](https://zaoqu-liu.github.io/SCENT/articles/algorithm.html) | Mathematical background and derivation |
| [Visualization](https://zaoqu-liu.github.io/SCENT/articles/visualization.html) | Publication-quality figures |
| [Benchmark](https://zaoqu-liu.github.io/SCENT/articles/benchmark.html) | Performance analysis |

## Citation

If you use SCENT in your research, please cite:

> Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of differentiation potency from a cell's transcriptome. *Nature Communications*. 2017;8:15599. doi:[10.1038/ncomms15599](https://doi.org/10.1038/ncomms15599)

```bibtex
@article{teschendorff2017single,
  title={Single-cell entropy for accurate estimation of differentiation 
         potency from a cell's transcriptome},
  author={Teschendorff, Andrew E and Enver, Tariq},
  journal={Nature Communications},
  volume={8},
  pages={15599},
  year={2017},
  doi={10.1038/ncomms15599}
}
```

## Authors

- **Zaoqu Liu** — Maintainer, C++ optimization ([GitHub](https://github.com/Zaoqu-Liu) · [ORCID](https://orcid.org/0000-0002-0452-742X))
- **Andrew E Teschendorff** — Original algorithm developer

## Links

- 📦 **R-Universe**: https://zaoqu-liu.r-universe.dev/SCENT
- 📖 **Documentation**: https://zaoqu-liu.github.io/SCENT/
- 🐛 **Issues**: https://github.com/Zaoqu-Liu/SCENT/issues
- 📄 **Original Paper**: [Nature Communications 2017](https://www.nature.com/articles/ncomms15599)
