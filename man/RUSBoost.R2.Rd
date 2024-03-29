% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{RUSBoost.R2}
\alias{RUSBoost.R2}
\title{RUSBoost variant of AdaBoost.R2}
\usage{
RUSBoost.R2(
  form,
  train,
  test,
  t_final = 100,
  power = 2,
  perc.U = 0.9,
  rel.thr = 0.9,
  coef = 1.5,
  ...
)
}
\arguments{
\item{form}{The model formula.}

\item{train}{A data.frame with the training data.}

\item{test}{A data.frame with the test data.}

\item{t_final}{The number of maximum boosting iterations. Default is 100.}

\item{power}{Type of loss function, e.g. linear (1), squared (2). Default is 2.}

\item{perc.U}{Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.}

\item{rel.thr}{Relevance threshold. Default is 0.9.}

\item{coef}{Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.}

\item{...}{Dots are passed to rpart}
}
\value{
Returns a vector with the predictions made by RUSBoost.R2.
}
\description{
RUSBoost variant of AdaBoost.R2
}
\examples{
data(Boston,package="MASS")

idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
form <- medv ~ .

train <- Boston[idx,]
test <- Boston[-idx,]

preds <- RUSBoost.R2(form,train,test)

}
