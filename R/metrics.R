# Mean reciprocal rank
mrr <- function(x) {
  x <- x[!is.na(x)]
  mean(1/x)
}
