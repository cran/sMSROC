\name{variance_probs-internal}
\alias{variance_probs
}
\title{Variance of the predictive model
}
\description{Estimation of the variance of the predictive model by bootstrap.
}
\usage{variance_probs(marker, outcome, status, observed.time, left, right, time,
               meth, data_type, grid, probs, ci.nboots, parallel, ncpus, all)
}
\arguments{
  \item{marker}{vector with the biomarker values.
}
  \item{outcome}{vector with the condition of the subjects as positive, negative or unknown at the considered time \code{time}.
}
  \item{status}{response vector with the outcome values. The highest one is assumed to stand for the subjects having the event under study.
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
  \item{probs}{vector containing the probabilities estimated through to the predictive model.
}
  \item{ci.nboots}{number of bootstrap samples.
}
  \item{parallel}{indicates whether parallel computing will be done or not.

}
  \item{ncpus}{number of CPUs to use if parallel computing is performed.
}
  \item{all}{indicates whether the probabilities from the predictive model should be considered or not.
}
}
\value{List with a single component:
\item{sd.probs }{vector containing the standard deviation of the probabilities of the predictive model.}
}


