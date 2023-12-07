\name{sMS_timeic-internal}
\alias{sMS_timeic}
\title{sMS estimator for prognostic biomarkers and interval censoring}
\description{Wrap function for computing the sMS estimator in prognosis scenarios under interval censorship.}
\usage{sMS_timeic(marker, left, right, outcome, time, meth, grid, probs, all)}
\arguments{
  \item{marker}{vector with the biomarker values.}
   \item{left}{vector containing the lower edges of the observed intervals.}
   \item{right}{vector with the upper edges of the observed intervals. The infinity is admissible as value (indicated as \strong{inf}).}
  \item{outcome}{ vector containing the condition of the individuals as positive, negative or censored at the time \code{time}.}
  \item{time}{point of time at which the sMSROC curve estimator will be computed.}
  \item{meth}{method for approximating the predictive model \eqn{P(D|X=x)}.
\itemize{
\item \bold{``E'',} allocates to each individual their own condition as positive or negative. Those whose condition is unknown at time \code{time} are dismissed.
\item \bold{``L'',} for proportional hazards regression models taking into account the observation intervals (see Details \code{\link{sMSROC}}).
\item \bold{``S'',} for proportional hazards regression models \bold{without} taking into account the observation intervals (see Details \code{\link{sMSROC}}).
}
}
  \item{grid}{grid size.}
  \item{probs}{vector with the probabilities from the predictive model when it is manually entered.}
\item{all}{parameter indicating whether all probabilities given by the predictive model should be considered (value \bold{``T''}) or just those corresponding to individuals whose condition as positive or negative is unknown (\bold{``F''}). The default value is (\bold{``T''}).
}
}

\details{This function gets the probabilities corresponding to the predictive model (first stage of the sMS ROC curve estimator). If they were not manually entered, the functions \code{pred.mod.emp} or \code{pred.mod.timeic} are called depending on the chosen \code{meth}. Then, it calls the function \code{computeROC} to compute the weighted empirical ROC curve estimator (second stage).}

\value{The returned value is a list with the following components:
\item{SE}{vector with the weighted empirical estimator of the sensitivity.}
\item{SP}{vector with the weighted empirical estimator of the specificity.}
\item{u}{vector containing the points between 0 and 1 at which the ROC curve estimator will be computed. Its size is determined by the \code{grid} parameter.}
\item{ROC}{ROC curve approximated at each point of the vector \code{u}.}
\item{auc}{area under the weighted empirical ROC curve estimator.}
\item{marker}{vector with the ordered biomarker values.}
\item{outcome}{vector with the condition of the individuals at time \code{time} as positive, negative or unknown.}
\item{probs}{vector with the probabilities of the predictive model corresponding to each biomarker value.}
}
\seealso{
\code{pred.mod.emp}, \code{pred.mod.binout},
\code{computeROC}, \code{sMS.binout} and \code{sMS.timerc}
}
