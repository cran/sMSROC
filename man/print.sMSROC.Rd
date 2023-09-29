\name{print.sMSROC}
\alias{print.sMSROC}
\title{Print sMSROC
}
\description{Prints the estimated AUC and the probabilistic model used to compute the predictive model.
}
\usage{
\method{print}{sMSROC}(x, ...)
}
\arguments{
  \item{x}{object of class \code{sMS} returned by the function \code{\link{sMSROC}}.
  }

  \item{...}{Ignored.}
}

\details{This function prints the estimated area under the ROC curve computed through the sMSROC estimator and the probabilistic model used to compute the predictive model.
}
\value{Printed output in the console containing the information described above.
}
\seealso{
\code{sMSROC}
}
\examples{
data(diabet)
roc <- sMSROC(marker=diabet$stab.glu, status=diabet$glyhb, conf.int="T")
print(roc)
}
