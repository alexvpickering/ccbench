% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_unique_drugs}
\alias{get_unique_drugs}
\title{Get unique drugs}
\usage{
get_unique_drugs(x)
}
\arguments{
\item{x}{Named vector where names are drug names.}
}
\value{
Unique drug names in order of appearance in \code{x}.
}
\description{
Used to compare results of queries where signatures are aggregated based on
the same drug or not.
}
