\name{pred_model_emp-internal}
\alias{pred_model_emp}
\title{Predictive model (naive estimation) }
\description{Naive estimation of the predictive model.}
\usage{pred_model_emp(marker, status)}
\arguments{
 \item{marker}{vector with the biomarker values.
 }
 \item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.
}
}
\details{This method for estimating the predictive model is used in both diagnosis and prognosis scenarios. It allocates individuals their own condition as positive or negative. Those with unknown condition are dismissed.
}

\value{The returned value is a list with two components:

\item{marker}{vector containing the ordered marker values.
}
\item{probs}{vector with the probabilities corresponding to each marker value estimated through the predictive model.
}
}
\seealso{
\code{sMS_binout}, \code{sMS_timerc}, \code{sMS_timeic} and \code{sMSROC}
}
