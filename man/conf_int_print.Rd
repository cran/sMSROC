\name{conf_int_print}
\alias{conf_int_print}
\title{AUC and confidence intervals
}
\description{Prints the AUC estimate value and its confidence intervals computed by the \code{\link{sMSROC}} function.
}
\usage{conf_int_print(sMS)
}
\arguments{
  \item{sMS}{object of class \code{sMS} returned by the function \code{\link{sMSROC}}.
}
}
\details{This function reads the AUC, lower and upper edges of its confidence intervals and the confidence level at which they were computed and prints this information in a single line.
}
\value{Printed string in the console containing the AUC, its confidence intervals and the confidence level at which they were computed.
}
\seealso{
\code{sMSROC}
}
\examples{
data(diabet)
roc <- sMSROC(marker=diabet$stab.glu, status=diabet$glyhb, conf.int="T")
conf_int_print(roc)
}
