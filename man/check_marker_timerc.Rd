\name{check_marker_timerc-internal}
\alias{check_marker_timerc}
\title{Check of prognosis scenarios under right censorship}
\description{
Checks the consistency of the parameters entered for prognosis scenarios under right censorship.
}
\usage{
check_marker_timerc(marker, status, observed.time, time, probs, sd.probs)
}

\arguments{
  \item{marker}{ vector with the biomarker values. It is a mandatory parameter.}
  \item{status}{ numeric response vector. Only two values will be taken into account. The highest one is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered. It is mandatory in prognosis scenarios and right censorship.}
  \item{observed.time}{ vector with the observed times for each subject. These values can be the event times or the censoring times. It is mandatory when dealing with time-dependent outcomes under right censorship.}
  \item{time}{ point of time at which the sMS ROC curve estimator will be computed. It is a mandatory parameter. The default value is 1.}
  \item{probs}{ vector containing the probabilities corresponding to the predictive model when it has been externally computed. Only values between [0,1] are admissible.}
  \item{sd.probs}{ vector with the standard deviations of the probabilities entered in \code{probs}. It is an optional parameter.}

}
\value{The function returns a list with the following components:
\item{marker}{ vector containing the biomarker values.}
\item{status}{ response vector.}
\item{observed.time}{ vector containing the observed time. Recall event/censoring time.}
\item{probs}{ vector with the probabilities corresponding to the predictive model.}
\item{sd.probs}{ vector containing the standard deviations of the predictive model if they have been manually entered.}
\item{outcome}{vector with the condition of the subjects at the time given in \code{time} as positive, negative or unknown.}
\item{controls}{ number of negative subjects.}
\item{cases}{ number of positive subjects.}
\item{misout}{ number of unknown subjects.}
\item{message}{ table containing the warning messages generated during the execution of the function.}
}

\seealso{
\code{check_marker_binout} and \code{check_marker_timeic}
}

