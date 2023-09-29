\name{auc_ci_empr-internal}
\alias{auc_ci_empr
}
\title{Confidence intervals for the AUC (empirical variance estimation)
}
\description{Computation of confidence intervals for the AUC by implementing the empirical
procedure for estimating the variance of the AUC, as described in \doi{10.1515/ijb-2019-0097}.
}
\usage{auc_ci_empr(SE, SP, auc, probs, controls, cases, ci.cl)
}
\arguments{
  \item{SE}{vector containing the values of the sensitivity returned from the \code{\link{sMSROC}} function.
}
  \item{SP}{vector containing the values of the specificity.
}
  \item{auc}{value with the AUC estimate.
}
  \item{probs}{vector containing the probabilities estimated through the predictive model.
}
  \item{controls}{number of negative individuals.}
  \item{cases}{number of positive individuals.}
  \item{ci.cl}{confidence level at which confidence intervals will be computed.
}
}
\value{List with two components:
  \item{ic.l}{lower edge of the confidence interval.}
  \item{ic.u}{upper edge of the confidence interval.}
}

