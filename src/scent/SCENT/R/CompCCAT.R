#' Correlation of Connectome And Transcriptome (CCAT)
#'
#' @description
#' Estimates signaling entropy via correlation between gene expression
#' and network connectivity.
#'
#' @param exp.m Gene expression matrix (genes x cells). Supports dense or sparse.
#' @param ppiA.m PPI network adjacency matrix (genes x genes).
#'
#' @return Numeric vector of CCAT scores (one per cell). Higher values indicate
#'   higher differentiation potency.
#'
#' @details
#' CCAT computes the Pearson correlation between gene expression and network
#' degree (number of interaction partners) for each cell. This approximation
#' is based on the observation that pluripotent cells express hub genes at
#' higher levels.
#'
#' The method automatically applies log2(x+1) transformation if the input
#' appears to be raw counts (max > 50).
#'
#' @author Andrew E Teschendorff, Zaoqu Liu
#' @references Teschendorff AE, Enver T. Nat Commun. 2017;8:15599.
#'
#' @examples
#' \dontrun{
#' data(net13Jun12.m)
#' exp <- matrix(rpois(5000 * 100, 5), nrow = 5000)
#' rownames(exp) <- head(rownames(net13Jun12.m), 5000)
#' ccat <- CompCCAT(exp, net13Jun12.m)
#' }
#'
#' @export
#' @importFrom stats cor var
#' @import Matrix
CompCCAT <- function(exp.m, ppiA.m) {
  
  # Input validation
  if (is.null(rownames(exp.m))) {
    stop("exp.m must have rownames (gene identifiers).")
  }
  if (is.null(rownames(ppiA.m))) {
    stop("ppiA.m must have rownames (gene identifiers).")
  }
  
  # Log-transform raw counts
  max.val <- max(exp.m, na.rm = TRUE)
  if (is.finite(max.val) && max.val > 50) {
    exp.m <- log2(exp.m + 1)
  }
  
  # Find common genes
  common.genes <- intersect(rownames(ppiA.m), rownames(exp.m))
  
  if (length(common.genes) < 5000) {
    stop("Insufficient gene overlap (<5000). Check gene identifiers.")
  }
  
  # Subset to common genes
  ppi.idx <- match(common.genes, rownames(ppiA.m))
  exp.idx <- match(common.genes, rownames(exp.m))
  
  ppi.sub <- ppiA.m[ppi.idx, ppi.idx]
  exp.sub <- exp.m[exp.idx, , drop = FALSE]
  
  # Compute network degree
  degree.v <- Matrix::rowSums(ppi.sub)
  
  # Check degree variance
  if (var(degree.v) < 1e-10) {
    warning("All genes have similar degree. CCAT may not be informative.")
    return(rep(NA_real_, ncol(exp.sub)))
  }
  
  # Convert sparse to dense for cor()
  if (inherits(exp.sub, "sparseMatrix")) {
    exp.sub <- as.matrix(exp.sub)
  }
  
  # Compute correlation between expression and degree
  ccat.v <- as.vector(cor(exp.sub, degree.v))
  
  return(ccat.v)
}
