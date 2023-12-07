\name{sMS_binout-internal}
\alias{sMS_binout}
\title{sMS estimator for diagnostic biomarkers}
\description{Wrap function for computing the sMS estimator in diagnosis scenarios.}
\usage{sMS_binout(marker, status, meth, grid, probs, all)}
\arguments{
  \item{marker}{ vector with the biomarker values.}
  \item{status}{ numeric response vector.}
  \item{meth}{ method for approximating the predictive model \eqn{P(D|X=x)}.
\itemize{
\item \bold{``E'',} allocates to each individual their own condition as positive or negative. Those whose condition is unknown at time \code{time} are dismissed.
\item \bold{``L'',} for \bold{L}inear logistic regression models (see details in \code{\link{sMSROC}}).

\item \bold{``S'',} for \bold{S}mooth models (see details \code{\link{sMSROC}}).}
}
\item{grid}{ grid size.
  }
\item{probs}{vector with the probabilities from the predictive model when it is manually entered.
}
\item{all}{parameter indicating whether all probabilities given by the predictive model should be considered (value \bold{``T''}) or just those corresponding to individuals whose condition as positive or negative is unknown (\bold{``F''}). The default value is (\bold{``T''}).
}
}

\details{ The function obtains the probabilities corresponding to the predictive model (first stage of the sMS ROC curve estimator). If they were not manually entered, the functions \code{pred.mod.emp} or \code{pred.mod.binout} are called depending on the chosen \code{meth}. Then, it calls the function \code{computeROC} to compute the weighted empirical ROC curve estimator (second stage).}

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
\code{pred_mod_emp}, \code{pred_mod_binout},
\code{computeROC}, \code{sMS_timerc} and \code{sMS_timeic}
}


