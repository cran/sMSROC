\name{auc_ci_nvar-internal}
\alias{auc_ci_nvar
}
\title{Confidence intervals for the AUC (theoretical variance estimation)
}
\description{Computation of confidence intervals for the AUC by implementing the theoretical procedure for estimating the variance of the AUC, as described in  \doi{10.1515/ijb-2019-0097}.
}
\usage{auc_ci_nvar(marker, outcome, status, observed.time, left, right, time,
                    meth, data_type, grid, probs, sd.probs, ci.cl, nboots,
                    SE, SP, auc, parallel, ncpus, all)
}
\arguments{
  \item{marker}{vector with the biomarker values.
}
  \item{outcome}{vector with the condition of the subjects as positive, negative or unknown at the considered time \code{time}.
}
  \item{status}{response vector.
}
  \item{observed.time}{vector with the observed times for each subject.
}
  \item{left}{vector with the lower edges of the observed intervals.
}
  \item{right}{vector with the upper edges of the observed intervals.
}
  \item{time}{point of time at which the sMS ROC curve estimator will be computed.
}
  \item{meth}{method for approximating the predictive model \eqn{P(D|X=x)}.
}
  \item{data_type}{scenario handled.
}
  \item{grid}{grid size.
}
  \item{probs}{vector containing the probabilities estimated through the predictive model.
}
  \item{sd.probs}{vector containing the standard deviation of the probabilities of the predictive model.
}
  \item{ci.cl}{confidence levet at which the confidence intervals will be computed.
}
  \item{nboots}{number of bootstrap samples.
}
  \item{SE}{vector containing the values of the sensitivity returned from \code{\link{sMSROC}} function.
}
  \item{SP}{vector containing the values of the specificity.
}
  \item{auc}{value with the AUC estimate.
}
  \item{parallel}{indicates whether parallel computing will be performed or not.
}
  \item{ncpus}{number of CPUs to use if parallel computing is performed.
}
  \item{all}{parameter indicating whether all probabilities given by the predictive model should be considered (value \bold{``T''}) or just those corresponding to individuals whose condition as positive or negative is unknown (\bold{``F''}). The default value is (\bold{``T''}).
}
}
\value{List with two components:
  \item{ic.l}{lower edge of the confidence interval.}
  \item{ic.u}{upper edge of the confidence interval.}
}

