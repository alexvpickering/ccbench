#' Mean reciprocal rank
#' @export
mrr <- function(x) {
  x <- x[!is.na(x)]
  mean(1/x)
}
