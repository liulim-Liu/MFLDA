% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvmflda.R
\name{cvMFLDA}
\alias{cvMFLDA}
\title{Cross-Validation for Multivariate Functional Linear Discriminant Analysis (MFLDA)}
\usage{
cvMFLDA(
  Xdata = Xdata,
  Y = Y,
  withCov = FALSE,
  plotIt = FALSE,
  metrics.choice = "Accuracy",
  Xtestdata = NULL,
  Ytest = NULL,
  isParallel = TRUE,
  ncores = NULL,
  gridMethod = "RandomSearch",
  AssignClassMethod = "Joint",
  nfolds = 5,
  ngrid = 8,
  standardize = TRUE,
  maxiteration = 20,
  weight = 0.5,
  thresh = 0.001
)
}
\arguments{
\item{Xdata}{A data frame containing the training data in long format, with the first three columns for `ID`, `Time`, and `Group`, followed by variable columns.}

\item{Y}{A vector representing the group labels for the training data.}

\item{plotIt}{A logical value indicating whether to generate discriminant plots. Default is FALSE.}

\item{metrics.choice}{A string specifying the metric to optimize. Options are "Accuracy" or "CombinedMetrics". Default is "Accuracy".}

\item{Xtestdata}{A data frame containing the test data. If NULL, training data is used for testing.}

\item{Ytest}{A vector representing the group labels for the test data. Required if Xtestdata is provided.}

\item{isParallel}{A logical value indicating whether to perform parallel processing. Default is TRUE.}

\item{ncores}{An integer specifying the number of cores to use for parallel processing. Default is NULL, which uses half the available cores.}

\item{nfolds}{An integer specifying the number of folds for cross-validation. Default is 5.}

\item{ngrid}{An integer specifying the number of tuning grid values. Default is 8.}

\item{standardize}{A logical value indicating whether to standardize the data to have mean zero and variance one for each time point and each variable. Default is TRUE.}

\item{maxiteration}{An integer specifying the maximum number of iterations for optimization. Default is 20.}

\item{thresh}{A numeric value indicating the convergence threshold. Default is 1e-03.}
}
\value{
A list containing:
  \item{CVOut}{A matrix of cross-validation results for different tuning parameter values.}
  \item{mfldaerror.test}{The estimated classification error for the test data.}
  \item{sidaerror.train}{The estimated classification error for the training data.}
  \item{hatalpha}{The estimated alpha coefficients from the MFLDA.}
  \item{PredictedClass}{The predicted class labels for the test data.}
  \item{var.selected}{The variables selected by the MFLDA.}
  \item{varname.selected}{The names of the selected variables.}
  \item{optTau}{The optimal tuning parameter values.}
  \item{gridValues}{The grid of tuning parameter values used.}
  \item{myDiscPlot}{A ggplot object of the discriminant plot if \code{plotIt} is TRUE; otherwise NULL.}
  \item{InputData}{The original input data used.}
}
\description{
This function performs cross-validation for Multivariate Functional Linear Discriminant Analysis (MFLDA).
It allows for optional parallel processing and can return various metrics, including accuracy or combined metrics.
}
\examples{
# Example usage of cvMFLDA
result <- cvMFLDA(Xdata, Y, plotIt=TRUE, metrics.choice="Accuracy",
                   Xtestdata=test_data, Ytest=test_labels,
                   nfolds=5, ngrid=8)

# Accessing results
print(result$mfldaerror.test)
print(result$myDiscPlot)

}
