#' Integrate Expression Data with PPI Network
#'
#' @description
#' Finds common genes between expression data and PPI network, then extracts
#' the maximally connected subnetwork.
#'
#' @param exp.m Gene expression matrix (genes x cells). Should be library-size
#'   normalized. Raw counts will be log2-transformed automatically.
#' @param ppiA.m PPI network adjacency matrix (genes x genes).
#'
#' @return List with elements:
#' \describe{
#'   \item{expMC}{Expression matrix for genes in the maximal component.}
#'   \item{adjMC}{Adjacency matrix of the maximal connected component.}
#' }
#'
#' @details
#' The function:
#' 1. Identifies genes present in both expression data and PPI network
#' 2. Extracts the induced subnetwork
#' 3. Finds the largest connected component
#' 4. Returns matched expression and adjacency matrices
#'
#' @author Andrew E Teschendorff, Zaoqu Liu
#' @references Teschendorff AE, Enver T. Nat Commun. 2017;8:15599.
#'
#' @examples
#' \dontrun{
#' data(net13Jun12.m)
#' exp <- matrix(rpois(5000 * 100, 5), nrow = 5000)
#' rownames(exp) <- head(rownames(net13Jun12.m), 5000)
#' integ <- DoIntegPPI(exp, net13Jun12.m)
#' }
#'
#' @export
#' @import Matrix
#' @import igraph
DoIntegPPI <- function(exp.m, ppiA.m) {
  
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
  
  # Extract subnetwork
  ppi.idx <- match(common.genes, rownames(ppiA.m))
  adj.m <- ppiA.m[ppi.idx, ppi.idx]
  
  # Convert to sparse for memory efficiency
  if (!inherits(adj.m, "dgCMatrix")) {
    adj.m <- Matrix::Matrix(adj.m, sparse = TRUE)
  }
  
  # Find maximal connected component
  gr <- igraph::graph_from_adjacency_matrix(adj.m, mode = "undirected")
  comp <- igraph::components(gr)
  
  if (comp$no == 0) {
    stop("Network has no connected components.")
  }
  
  # Find largest component
  largest.id <- which.max(comp$csize)
  mc.idx <- which(comp$membership == largest.id)
  
  if (length(mc.idx) < 1000) {
    warning("Largest connected component has <1000 genes.")
  }
  
  # Extract maximal component
  adjMC.m <- adj.m[mc.idx, mc.idx]
  
  # Match expression matrix
  exp.idx <- match(rownames(adjMC.m), rownames(exp.m))
  expMC.m <- exp.m[exp.idx, , drop = FALSE]
  
  list(expMC = expMC.m, adjMC = adjMC.m)
}
