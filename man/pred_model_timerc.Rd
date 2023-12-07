\name{pred_model_timerc-internal}
\alias{pred_model_timerc}
\title{Predictive model in prognosis scenarios (II)}
\description{Estimation of the predictive model in prognosis scenarios under right censorship.}
\usage{pred_model_timerc(marker, status, observed.time, outcome, time, meth)}
\arguments{
\item{marker}{vector with the biomarker values.}
\item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.
}
\item{observed.time}{ vector with the observed times for each subject. Notice that these values may be the event times or the censoring times.}
\item{outcome}{vector with the status of the subjects as positive, negative or censored (unknown) at the considered time \code{time}.}
\item{time}{point of time at which the sMS ROC curve estimator will be computed.}
\item{meth}{method for approximating the predictive model \eqn{P(D|X=x)}. The options are:
\itemize{
\item \bold{``L'',} for proportional hazards regression models (see Details).
\item \bold{``S'',} for smooth models (see Details).
}}
}

\details{
\itemize{
\item If \code{meth} = \bold{``L''}, the event times are assumed to come from a Cox proportional hazards regression model:
 \deqn{ P (T \leq t \;|\; X=x) = 1 - \exp \{ - \Delta_0(t) \cdot \exp \{ \beta_0 + \beta_1 \cdot \log(x)\}\},} where \eqn{\Delta_0(\cdot)} is the baseline hazard function and \eqn{\beta_0, \beta_1 \in {\cal R}}.
\item If \code{meth} = \bold{``S''}, the approximation is done by \deqn{P (T\leq t \;|\; X=x) = 1 - \exp \{ - \Delta_0(t) \cdot \exp \{ s(x)\}\} } being \eqn{s(\cdot)} the smooth function (penalized splines, \doi{10.1111/1467-9868.00125}).}
}
\value{The returned value is a list with three components:
\item{marker}{vector containing the ordered marker values.}
\item{probs}{vector with the probabilities corresponding to each marker value estimated through the predictive model.}
\item{outcome}{vector with the status of the subjects as positive, negative or censored at the considered time \code{time}.}
}

\seealso{
\code{sMS_timerc} and \code{sMSROC}
}

