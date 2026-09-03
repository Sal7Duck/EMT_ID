#' Protein-Protein Interaction Network
#'
#' @description
#' Protein-protein interaction network adjacency matrix compiled from
#' multiple public databases. Contains 8434 genes.
#'
#' @format A symmetric binary matrix (8434 x 8434):
#' \describe{
#'   \item{Rows/Columns}{Gene Entrez IDs}
#'   \item{Values}{1 = interaction, 0 = no interaction}
#' }
#'
#' @source Compiled from public PPI databases (June 2012).
#' @references Teschendorff AE, Enver T. Nat Commun. 2017;8:15599.
#'
#' @examples
#' data(net13Jun12.m)
#' dim(net13Jun12.m)
#' sum(net13Jun12.m) / 2  # Number of interactions
#'
"net13Jun12.m"
