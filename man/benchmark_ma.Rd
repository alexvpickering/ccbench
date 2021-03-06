% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{benchmark_ma}
\alias{benchmark_ma}
\title{Test meta-analysis against contrasts for a drug in CMAP.}
\usage{
benchmark_ma(
  dprimes,
  anals,
  drug_info = NULL,
  suffix = "",
  data_dir = getwd(),
  pvals = NULL,
  pval_breaks = c(seq(0.01, 0.3, 0.01), seq(0.3, 1, 0.05))
)
}
\description{
Test meta-analysis against contrasts for a drug in CMAP.
}
