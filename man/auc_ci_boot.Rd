\name{auc_ci_boot-internal}
\alias{auc_ci_boot
}
\title{Confidence intervals for the AUC (bootstrap)
}
\description{Computation of confidence intervals for the AUC based on Bootstrap Percentile.
}
\usage{auc_ci_boot(marker, outcome, status, observed.time, left, right, time,
                    data_type, meth, grid, probs, ci.cl, ci.nboots, parallel,
                    ncpus, all)
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
  \item{data_type}{scenario handled.
}
  \item{meth}{method for approximating the predictive model \eqn{P(D|X=x)}.
}
  \item{grid}{grid size.
}
  \item{probs}{vector containing the probabilities estimated through the predictive model.
}
  \item{ci.cl}{confidence level at which the confidence intervals will be computed.
}
  \item{ci.nboots}{number of bootstrap samples.
}
  \item{parallel}{indicates whether parallel computing will be performed or not.

}
  \item{ncpus}{number of CPUs to use if parallel computing is performed.
}
  \item{all}{indicates whether the probabilities from the predictive model will be considered for all individuals, or only for those whose outcome value (condition) is unknown.
}
}
\value{List with two components:
  \item{ic.l}{ lower edge of the confidence interval.}
  \item{ic.u}{ upper edge of the confidence interval.}

}
