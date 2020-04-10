#' Get unique drugs
#'
#' Used to compare results of queries where signatures are aggregated based on
#' the same drug or not.
#'
#' @param x Named vector where names are drug names.
#'
#' @return Unique drug names in order of appearance in \code{x}.
#' @export
#'
get_unique_drugs <- function(x) {
  drugs <- gsub('^([^_]+)_.+?$', '\\1', names(x))
  unique(drugs)
}


#' Extract Adjusted P-Values
#'
#' Used to benchmark p-value cutoffs for query signatures.
#'
#' @param es result of call to \code{es_meta}.
#' @param anals result of call to \code{diff_expr)}.
#'
#' @return List with slots:
#'   \item{meta}{Named numeric vector with overall false discovery rates for all gene from meta-analysis.}
#'   \item{contrasts}{List of named numeric vectors (one per contrast) with limma adjusted p-values for all genes.}
#' @export
#' @seealso \code{\link[crossmeta]{es_meta}}, \code{\link[crossmeta]{diff_expr}}.
#'
get_pvals <- function(es, anals) {

  pvals <- list()
  tts <- lapply(anals, `[[`, 'top_tables')
  tts <- unname(tts)
  tts <- unlist(tts, recursive = FALSE)

  for (i in seq_along(es)) {
    namei <- names(es)[i]
    esi <- es[[i]]
    meta <- esi$filt$fdr
    names(meta) <- row.names(esi$filt)

    contrasts <- list()
    dp_cols <- grep("^dp", colnames(esi$raw), value = TRUE)

    for (col in dp_cols){
      con <-  gsub('^dprime.', '', col)
      pvalsi <- tts[[con]]$adj.P.Val
      names(pvalsi) <- row.names(tts[[con]])
      contrasts[[con]] <- pvalsi[!is.na(pvalsi)]
    }
    pvals[[namei]] <- list(meta=meta, contrasts=contrasts)
  }
  return(pvals)
}


#' Get ranks for query result
#'
#' @param query_res Lists of query results with names corresponding to query signature. Each list
#'   contains a sorted numeric vector of pearson correlation with names corresponding to queried signatures.
#'
#' @return Numeric vector of ranks with the same names as \code{query_res}. NA indicates no query
#'    took place (e.g. if no genes met p-value threshold for query).
#' @export
#'
#' @examples
get_ranks <- function(query_res) {

  query_names <- names(query_res)
  ranks <- c()
  for (i in seq_along(query_res)) {
    x <- query_res[[i]]
    x <- get_unique_drugs(x)
    query <- query_names[i]

    ranks[query] <- if(length(x)) which(x == query) else NA
  }
  names(ranks) <- query_names
  return(ranks)
}
