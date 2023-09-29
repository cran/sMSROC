\name{check_marke_timeic-internal}
\alias{check_marker_timeic}
\title{Check of prognosis scenarios under interval censorship
}
\description{
Checks the consistency of the parameters entered for prognosis scenarios under interval censorship.
}
\usage{
check_marker_timeic(marker, left, right, time, probs, sd.probs)
}
\arguments{
  \item{marker}{ vector with the biomarker values. It is a mandatory parameter.}
  \item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios and interval censorship and ignored in other situations.}
  \item{right}{vector with the upper edges of the observed intervals. It is mandatory as well in prognosis scenarios and interval censorship and ignored in other situations.}
  \item{time}{ point of time at which the sMS ROC curve estimator will be computed. It is a mandatory parameter. The default value is 1.}
  \item{probs}{ vector containing the probabilities corresponding to the predictive model when it has been externally computed. Only values between [0,1] are admissible.}
  \item{sd.probs}{ vector with the standard deviations of the probabilities entered in \code{probs}. It is an optional parameter.}
}
\value{The function returns a list with the following components:
\item{marker}{ vector containing the biomarker values.}
\item{left}{ vector containing the lower edges of the observed intervals.}
\item{right}{ vector containing the upper edges of the observed intervals.}
\item{probs}{ vector with the probabilities corresponding to the predictive model.}
\item{sd.probs}{ vector containing the standard deviations of the predictive model if they have been manually entered.}
\item{outcome}{vector with the condition of the subjects at the time \code{time}, as positive, negative or unknown.}
\item{controls}{ number of negative subjects.}
\item{cases}{ number of positive subjects.}
\item{misout}{ number of subjects whose condition is unknown.}
\item{message}{ table containing the warning messages generated during the execution of the function.}
}

\seealso{
\code{check_marker_binout} and \code{check_marker_timerc}
}
