#' Query using p-values to choose number of genes
#'
#' @param dprimes List with names of signatures to query with. Each item of list is named a numeric vector of
#'   effect size values with names corresponding to the gene.
#' @param drug_info Matrix of differential expression values to query against.
#'   Row names are genes and column names are signatures.
#'
#' @param pvals Same as \code{dprimes} but p-values or FDRs in place of effect size values.
#' @param breaks Breaks to use for determing \code{ngenes} argument to \code{query_drugs}.
#'
#' @return Named list with items:
#'  \item{resl}{For each value of \code{breaks}, a list of pearson correlations between each signature in
#'  \code{dprimes} and \code{drug_info}.}
#'  \item{ngenes}{For each value of \code{breaks}, a numeric vector indicating the number of genes used for the query.}
#'
#'
#'  @seealso \code{\link[ccmap]{query_drugs}}
#'
#' @export
#'
query_pval_ngenes <- function(dprimes, drug_info, pvals, breaks = c(seq(0.01, 0.3, 0.01), seq(0.3, 1, 0.05))) {
  resl <- list()
  ngenes <- list()

  for (p in breaks) {
    cat('Working on pval:', p, '\n')
    pc <- as.character(p)
    ng <- sapply(pvals, function(con) sum(con <= p))

    res <- lapply(seq_along(dprimes), function(i) {
      con <- dprimes[[i]]
      ccmap::query_drugs(con, drug_info, ngenes = ng[i])
    })

    names(res) <- names(dprimes)
    resl[[pc]] <- res
    ngenes[[pc]] <- ng
  }

  return(list(resl = resl, ngenes = ngenes))
}


#' Query using fixed number of genes
#'
#' @inheritParams query_pval_ngenes
#' @inherit query_pval_ngenes return seealso
#'
#' @export
#'
query_fixed_ngenes <- function(dprimes, drug_info, breaks = c(seq(50, 950, 50), 1001)) {
  resl <- list()
  ngenes <- list()

  for (n in breaks) {
    cat('Working on ngenes:', n, '\n')
    nc <- as.character(n)
    ng <- sapply(dprimes, function(x) n)

    res <- lapply(seq_along(dprimes), function(i) {
      con <- dprimes[[i]]
      ccmap::query_drugs(con, drug_info, ngenes = ng[i])
    })

    names(res) <- names(dprimes)
    resl[[nc]] <- res
    ngenes[[nc]] <- ng
  }

  return(list(resl = resl, ngenes = ngenes))
}

