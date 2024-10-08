% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratUtilities.R
\name{fitGMM}
\alias{fitGMM}
\title{Rescale data using GMM}
\usage{
fitGMM(
  seu,
  method = c("Rmixmod", "mclust"),
  model = mixmodGaussianModel(family = "general", listModels = "Gaussian_p_Lk_Ck",
    free.proportions = FALSE, equal.proportions = TRUE),
  max_zero_percentage = 0.1,
  ...
)
}
\arguments{
\item{seu}{A Seurat object}

\item{method}{The package used to do the GMM fit. Available tools are
'Rmixmod' and 'mclust'.}

\item{model}{If method is Rmixmod, the model will be used to fit the data.
see \link[Rmixmod]{mixmodGaussianModel}.}

\item{max_zero_percentage}{The cutoff of zero percentage. If the percentage
is greater than max_zero_percentage, the zeros will be removed by random
pick to meet the max percentage. This is used to make sure there is enough
variance to fit the model.}

\item{...}{Not use.}
}
\value{
A Seurat object with a new assay named as 'GMM'. The re-scaled values
saved in layer 'data'.
}
\description{
Re-scale data using multiclass Gaussian Mixture Model. If data is lack of
variance, the function will trying to do log normal distribution re-scale.
Re-scaled values will be probabilities between 0-1.
}
