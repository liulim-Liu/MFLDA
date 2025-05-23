% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spline.R
\name{spline.prediction}
\alias{spline.prediction}
\title{Spline-Based Estimation from Dataset}
\usage{
spline.prediction(olddat)
}
\arguments{
\item{olddat}{A data frame containing the dataset with participant observations. It must include columns for `id`, `group`, `time`, and additional variables to be predicted.}
}
\value{
A list containing:
  \item{pred.df}{A data frame with the predicted values for each variable across time points, along with participant IDs and group classifications.}
}
\description{
This function takes a dataset of participant observations and uses spline basis functions to generate fitted values for each variable across unique time points. It constructs a linear model for each variable and predicts values using B-splines.
}
\details{
The function internally identifies unique time points from the dataset and generates B-spline basis matrices for these time points.
It uses linear models to predict values for each variable, storing the fitted values in a list.
The function also assigns group classifications based on participant IDs if the original data has irregular time points.
}
\examples{
# Assume `my_data` is a data frame with the necessary structure
prediction_results <- spline.prediction(my_data)

# View the predicted values
head(prediction_results$pred.df)

}
