\name{evol_auc}
\alias{evol_auc}
\title{Evolution of the AUCs
}
\description{
Plots, in prognosis scenarios, the areas under the ROC curves computed by the sMSROC estimator for a sequence of times.
}
\usage{evol_auc(marker, status, observed.time, left, right,
         time = 1, meth = c("L", "S", "E"), grid = 500)
}
\arguments{
  \item{marker}{ vector with the biomarker values.
}
  \item{status}{  numeric response vector. Only two values will be taken into account. The highest one is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered. It is mandatory in prognosis scenarios and right censorship.
}
  \item{observed.time}{ vector with the observed times for each subject, for prognosis scenarios under right censorship. Notice that these values may be the event times or the censoring times.
}
 \item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations.
}
  \item{right}{vector with the upper edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations. The infinity is admissible as value (indicated as \bold{inf}).
}
  \item{time}{vector of times at which the sMS ROC curve estimator will be computed. The default value is 1.}

  \item{meth}{method for approximating the predictive model \eqn{P(D|X=x)}. There are several options available:

\itemize{
\item \bold{``E'',} allocates to each individual their own status value as having the event of interest or not. Those with missing status values or censored at a fixed point of time \bold{t} are dismissed.

\item \bold{``L'',} for proportional hazards regression models.

\item \bold{``S'',} for smooth models.
}
}
\item{grid}{grid size for computing the AUC. Default value 500.
}
}
\details{
This function calls the {\code{\link{sMSROC}}} function at each of the times indicated in the vector \code{time}, and the AUC is computed according to the parameters indicated.
}
\value{A list with the following components:
\item{evol.auc}{object of class \code{\link{ggplot}}. A graphic line plotting the AUCs at the considered times.}
\item{time}{vector with the ordered values of the \code{time} entered as parameter.}
\item{auc}{vector with the values of the AUCs computed at the times indicated at the \code{time} parameter.}
}

\seealso{
\code{sMSROC}
}
\examples{
# Example of the use of the evol.AUC function
data(ktfs)
DT = ktfs
aucs <- evol_auc(marker = DT$score,
                 status = DT$failure,
                 observed.time = DT$time,
                 time = seq(2:3),
                 meth = "E")
aucs$evol.auc
}
