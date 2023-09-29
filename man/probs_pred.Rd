\name{probs_pred}
\alias{probs_pred}
\title{Plot of the predictive model}
\description{This function plots the predicted probabilities for each marker value computed through the predictive model together, with 95\% pointwise confidence intervals.
}
\usage{probs_pred(sMS, var, nboots, parallel, ncpus)
}

\arguments{
  \item{sMS}{object of class \code{sMS} returned from function \code{\link{sMSROC}}.
}
  \item{var}{parameter indicating whether 95\% pointwise confidence intervals for the predictive model will be plotted (value \bold{"T"}) or not (value \bold{"F"}). The default value is "F".
}
  \item{nboots}{number of bootstrap samples to be generated for computing the pointwise confidence intervals. The default value is 500.
}
  \item{parallel}{parameter indicating whether parallel computing will be performed (value \bold{"T"}) or not (value \bold{"F"}). The default is ``F''.
}
  \item{ncpus}{number of CPUS to be used in the case of carrying out parallel computing. The default value is 1 and the maximum is 2.
}
}
\details{The function plots the probability function estimation of the predictive model versus the biomaker. It also computes and plots 95\% pointwise confidence intervals on the same graphic when the \code{var} parameter is set to "T".

The variance of the probability estimates, obtained by the predictive model, is computed via bootstrap with \code{nboots} samples.
}
\value{A list with these components:
\item{plot }{object of class \code{\link{ggplot}} (graphical output).}
\item{thres }{ordered biomarker values (x-axis coordinates).}
\item{probs }{predicted probabilities (y-axis coordinates).}
\item{sd.probs }{estimates of the standard deviation of the predicted probabilities.}
}
\seealso{\code{pred_model_binout}, \code{pred_model_timerc} and \code{pred_model_timeic}
}
\examples{
data(ktfs)
DT <- ktfs
roc <- sMSROC(marker = DT$score,
              status = DT$failure,
              observed.time = DT$time,
              time = 5,
              meth = "S")
probs <- probs_pred(sMS = roc)
probs$plot
}

