% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{normalize.loess}
\alias{normalize.loess}
\title{loess normalization}
\usage{
normalize.loess(
  mat,
  subset = sample(seq.int(nrow(mat)), min(c(5000, nrow(mat)))),
  span = 2/3,
  family.loess = "symmetric"
)
}
\arguments{
\item{mat}{a matrix with columns containing the values of the chips
to normalize.}

\item{subset}{a subset of the data to fit a loess to.}

\item{span, family.loess}{parameter to be passed the function
\link[stats]{loess}.}
}
\value{
the mean matrix to be subtracted.
}
\description{
loess normalization
}
