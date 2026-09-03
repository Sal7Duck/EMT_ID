#' Compute Signaling Entropy Rate (SR)
#'
#' @description
#' Computes single-cell differentiation potency using the signaling entropy
#' rate on a protein interaction network.
#'
#' @param integ.l Output from \code{\link{DoIntegPPI}}, containing expMC and adjMC.
#' @param local Logical. If TRUE, returns normalized local entropies. Default FALSE.
#'
#' @return List with elements:
#' \describe{
#'   \item{SR}{Signaling entropy rate (0-1). Higher = higher potency.}
#'   \item{inv}{Stationary distribution matrix (genes x cells).}
#'   \item{locS}{Local entropy matrix (genes x cells).}
#'   \item{nlocS}{Normalized local entropies, or NULL if local = FALSE.}
#' }
#'
#' @details
#' The signaling entropy rate is defined as:
#' \deqn{SR = \sum_j \pi_j S_j / SR_{max}}
#'
#' where \eqn{\pi_j} is the stationary distribution and \eqn{S_j} is the local
#' entropy at gene j. \eqn{SR_{max}} is the maximum entropy rate of the network,
#' computed from the largest eigenvalue of the adjacency matrix.
#'
#' @author Andrew E Teschendorff, Zaoqu Liu
#' @references Teschendorff AE, Enver T. Nat Commun. 2017;8:15599.
#'
#' @examples
#' \dontrun{
#' data(net13Jun12.m)
#' exp <- matrix(rpois(5000 * 50, 5), nrow = 5000)
#' rownames(exp) <- head(rownames(net13Jun12.m), 5000)
#' integ <- DoIntegPPI(exp, net13Jun12.m)
#' sr <- CompSRana(integ)
#' }
#'
#' @export
#' @useDynLib SCENT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#' @import igraph
CompSRana <- function(integ.l, local = FALSE) {
  
  # Input validation
  if (!is.list(integ.l) || !all(c("expMC", "adjMC") %in% names(integ.l))) {
    stop("integ.l must be output from DoIntegPPI (list with expMC and adjMC).")
  }
  
  exp.m <- integ.l$expMC
  adj.m <- integ.l$adjMC
  
  # Convert sparse to dense for C++
  if (inherits(exp.m, "sparseMatrix")) {
    exp.m <- as.matrix(exp.m)
  }
  if (inherits(adj.m, "sparseMatrix")) {
    adj.m <- as.matrix(adj.m)
  }
  
  # Compute maximum entropy rate
  maxSR <- .CompMaxSR(adj.m)
  
  if (!is.finite(maxSR) || maxSR <= 0) {
    stop("Invalid maxSR. Check if the network is connected.")
  }
  
  # Compute SR via C++
  result <- CompSRana_cpp(exp.m, adj.m, maxSR)
  
  # Normalized local entropy
  if (local) {
    degrees <- rowSums(adj.m > 0)
    # For degree <= 1, max entropy is 0, set nlocS to 0
    log.degrees <- log(degrees)
    log.degrees[degrees <= 1] <- Inf  # Will make nlocS = 0 for these
    nlocS.m <- result$locS / log.degrees
    nlocS.m[!is.finite(nlocS.m)] <- 0
  } else {
    nlocS.m <- NULL
  }
  
  list(
    SR = result$SR,
    inv = result$inv,
    locS = result$locS,
    nlocS = nlocS.m
  )
}


#' Compute Maximum Entropy Rate
#' @param adj.m Adjacency matrix.
#' @return Maximum entropy rate (log of largest eigenvalue).
#' @keywords internal
.CompMaxSR <- function(adj.m) {
  
  # Use ARPACK to find largest eigenvalue
  fa <- function(x, extra = NULL) {
    as.vector(adj.m %*% x)
  }
  
  ap.o <- igraph::arpack(
    fa,
    options = list(n = nrow(adj.m), nev = 1, which = "LM"),
    sym = TRUE
  )
  
  lambda <- ap.o$values[1]
  
  if (lambda <= 0) {
    stop("Largest eigenvalue is non-positive. Network may be disconnected.")
  }
  
  log(lambda)
}
