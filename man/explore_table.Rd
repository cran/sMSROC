\name{explore_table}
\alias{explore_table}
\title{Exploratory data analysis}
\description{This function provides descriptive statistics for the pooled sample and the samples of positive, negative individuals and those whose condition is unknown.}
\usage{explore_table(marker, status, observed.time, left, right, time, d, ...)}
\arguments{
  \item{marker}{vector with the biomarker values. It is a mandatory parameter.}

  \item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.}

  \item{observed.time}{vector with the observed times for each subject, for prognosis scenarios under right censorship. Notice that these values may be the event times or the censoring times.}

 \item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations.}

  \item{right}{vector with the upper edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations. The infinity is admissible as value (indicated as \strong{inf}).}

  \item{time}{point of time at which the sMS ROC curve estimator will be computed. The default value is 1.}
  \item{d}{number of decimal figures to which the results will be rounded.}
  \item{...}{additional parameters of the \code{flextable} function which allow to customize the output table.}}

\details{The function computes the following descriptive statistics for the pooled sample and the samples of the different groups of individuals: minimum, maximun, mean, variance, standard deviation, and first, second and third quartiles.}

\value{The ouput is a list with two components:
\item{summary}{ matrix whose columns are the statistics described above and the rows show the corresponding results for each sample.}
\item{table}{object of class \code{\link{flextable}} that represent the matrix in \code{summary} in a customizable table.}

}

\seealso{
\code{explore_plot}
}

\examples{
data(diabet)
explore_table(marker=diabet$stab.glu, status=diabet$diab)
}

