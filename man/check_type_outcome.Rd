\name{check_type_outcome-internal}
\alias{check_type_outcome}
\title{Check the type of scenario (diagnosis/prognosis)
}
\description{
Determines the type of scenario handled: diagnosis or prognosis, under right or interval censorship, according to the parameters entered by the user.
}
\usage{check_type_outcome(status, observed.time, left, right)}

\arguments{
  \item{status}{numeric response vector. Only two values will be taken into account. The highest one is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.}
  \item{observed.time}{vector with the observed times for each subject, when dealing with time-dependent outcomes under right censorship. These values can be the event times or the censoring times.}
  \item{left}{vector containing the lower edges of the observed intervals. It is mandatory when dealing with prognosis scenarios and interval censorship, and will be ignored in other situations.}
  \item{right}{vector with the upper edges of the observed intervals. It is mandatory as well in prognosis scenarios and interval censorship and ignored in other situations.}
}
\details{If both the vectors status and observed time are indicated the funtion assumes a prognosis scenario and right censorship. When only the vector status is entered, a diagnosis scenario is set up. If none of these parameters are indicated but the left and right ones, a prognosis scenario and interval censorship is assumed. Any other case, the function is not able to determine the type of scenario.

}
\value{
A list with a single component:
\item{type.outcome}{ string of length 6 with the following values:
\itemize{
  \item{\strong{"binout",}}{ in the case of diagnosis scenarios.}
  \item{\strong{"timerc",}}{ for prognosis scenarios and right censorship.}
  \item{\strong{"timeic",}}{ for prognosis scenarios and interval censorship.}
  \item{\strong{"unknow",}}{ if it is not possible to determine the type of scenario.}
  }
}
}
