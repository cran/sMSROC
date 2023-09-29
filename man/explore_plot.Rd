\name{explore_plot}
\alias{explore_plot}
\title{Graphical exploratory data analysis}
\description{Plots the kernel density estimations of the biomarker distributions on positive and negative individuals.}
\usage{explore_plot(marker, status, observed.time, left, right, time)}
\arguments{
  \item{marker}{ vector with the biomarker values.}

  \item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.}

  \item{observed.time}{ vector with the observed times for each subject, for prognosis scenarios under right censorship. Notice that these values may be the event times or the censoring times.}

 \item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations.}

  \item{right}{vector with the upper edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations. The infinity is admissible as value (indicated as \strong{inf}).}

  \item{time}{point of time at which the sMS ROC curve estimator will be computed. The default value is 1.}
}

\value{The ouput is a list with three components:
\item{plot}{object of class \code{ggplot}.}
\item{neg}{vector with the biomarker values on negative individuals.}
\item{pos}{vector with the marker values on positive individuals.}
}

\seealso{
\code{explore_table}
}
\examples{
data(diabet)
explore_plot(marker=diabet$stab.glu, status=diabet$diab)
}



