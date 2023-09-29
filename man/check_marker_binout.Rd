\name{check_marker_binout-internal}
\alias{check_marker_binout}
\title{Checks of diagnosis scenarios}
\description{
Checks the consistency of the parameters entered for diagnosis scenarios.}
\usage{
check_marker_binout(marker, status, probs, sd.probs)
}
\arguments{
  \item{marker}{ vector with the biomarker values. It is a mandatory parameter.
}
  \item{status}{ numeric response vector. Only two values will be taken into account. The highest one is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered. It is a mandatory parameter in diagnosis scenarios.
 }
  \item{probs}{ vector containing the probabilities corresponding to the predictive model when it has been externally computed. Obviosly, only values between [0,1] are admissible.}
  \item{sd.probs}{ vector with the standard deviations of the probabilities entered in \code{probs}. It is an optional parameter.}
}

\value{The ouput is a list with the following components:
\item{marker}{vector containing the biomarker values.}
\item{outcome}{vector with the condition of the subjects as positive or negative.}
\item{probs}{vector with the probabilities corresponding to the predictive model.}
\item{sd.probs}{ vector containing the standard deviations of the predictive model if they have been manually entered.}
\item{controls}{ number of negative subjects.}
\item{cases}{ number of negative subjects.}
\item{misout}{ number of subjects whose outcome value is not known.}
\item{message}{ table containing the warning messages generated during the execution of the function.}
}

\seealso{
\code{check_marker_timerc} and \code{check_marker_timeic}
}
