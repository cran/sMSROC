\name{compute_ROC-internal}
\alias{compute_ROC}
\title{Weighted empirical ROC curve estimator
}
\description{Computes the weighted empirical ROC curve estimator associated to the input biomarker.
}
\usage{compute_ROC(marker, probs, grid)}
\arguments{
  \item{marker}{vector with the biomarker values.}
  \item{probs}{vector containing the probabilities corresponding to the predictive model.}
  \item{grid}{grid size.}
}
\details{ This function computes the weighted empirical estimators for the sensitivity (SE) and specificity (SP) using as weights the probabilities given by the predictive model. Then, the ROC curve is approximated through linear interpolation of 1 - SP and SE and computed at a partition of the \eqn{[0,1]} interval of size \code{grid}.
}
\value{ The returned value is a list with the following components:
\item{SE}{ vector with the weighted empirical estimator of the sensitivity.}
\item{SP}{ vector with the weighted empirical estimator of the specificity.}
\item{u}{ vector containing the points between 0 and 1 at which the ROC curve estimator will be computed. Its size is determined by the \code{grid} parameter.}
\item{ROC}{ ROC curve approximated at each point of the vector \code{u}.}
\item{auc}{ area under the weighted empirical ROC curve estimator.}
\item{marker}{ vector with the ordered biomarker values.}
\item{probs}{ vector with the probabilities of the predictive model corresponding to each biomarker value.}
}
\seealso{
\code{sMSbinout}, \code{sMStimerc} and \code{sMStimeic}
}
