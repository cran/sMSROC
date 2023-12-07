\name{sMS_timerc}
\alias{sMS_timerc}
\title{sMS estimator for prognostic biomarkers and right censorship}
\description{Wrap function for computing the sMS estimator in prognosis scenarios under right censorship.}
\usage{sMS_timerc(marker, status, observed.time, outcome, time,
           meth, grid, probs, all)}
\arguments{
  \item{marker}{vector with the biomarker values.}
  \item{status}{numeric response vector.}
  \item{observed.time}{vector with the observed times. These values may be the event times or the censoring times.}
  \item{outcome}{vector containing the condition of the individuals as positive, negative or censored (unknown) at the time \code{time}.}
  \item{time}{point of time at which the sMSROC curve estimator will be computed.}
  \item{meth}{ method for approximating the predictive model \eqn{P(D|X=x)}.
\itemize{
\item \bold{``E'',} allocates to each individual their own condition as positive or negative. Those whose condition is unknown at time \code{time} are dismissed.
\item \bold{``L'',} for \bold{L}inear proportional hazards regression models (see details in \code{\link{sMSROC}}).
\item \bold{``S'',} for \bold{S}mooth models (see details in \code{\link{sMSROC}}).
}
}
\item{grid}{grid size.
}
\item{probs}{vector with the probabilities from the predictive model when it is manually entered.
}
\item{all}{parameter indicating whether all probabilities given by the predictive model should be considered (value \bold{``T''}) or just those corresponding to individuals whose condition as positive or negative is unknown (\bold{``F''}). The default value is (\bold{``T''}).
}
}

\details{ This function gets the probabilities corresponding to the predictive model (first stage of the sMS ROC curve estimator). If they were not manually entered, the functions \code{pred.mod.emp} or \code{pred.mod.timerc} are called depending on the chosen \strong{meth}. Then, it calls the function \code{computeROC} to compute the weighted empirical ROC curve estimator (second stage).}

\value{ The returned value is a list with the following components:
\item{SE}{ vector with the weighted empirical estimator of the sensitivity.}
\item{SP}{ vector with the weighted empirical estimator of the specificity.}
\item{u}{ vector containing the points between 0 and 1 at which the ROC curve estimator will be computed. Its size is determined by the \code{grid} parameter.}
\item{ROC}{ ROC curve approximated at each point of the vector \code{u}.}
\item{auc}{ area under the weighted empirical ROC curve estimator.}
\item{marker}{ vector with the ordered biomarker values.}
\item{outcome}{ vector with the condition of the individuals at time \code{time} as positive, negative or cenrored (unknown).}
\item{probs}{ vector with the probabilities of the predictive model corresponding to each biomarker value.}
}
\seealso{
\code{pred.mod.emp}, \code{pred.mod.binout},
\code{computeROC}, \code{sMS.timerc}, \code{sMS.timeic}
}
