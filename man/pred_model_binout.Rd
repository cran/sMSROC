\name{pred_model_binout-internal}
\alias{pred_model_binout}
\title{Predictive model estimation in diagnosis scenarios}
\description{Estimation of the predictive models in diagnosis scenarios.
}
\usage{pred_model_binout(marker, status, meth)}
\arguments{
 \item{marker}{vector with the biomarker values.
 }
  \item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest value, for those who do not. Any other value will not be considered.
}

\item{meth}{ method for approximating the predictive model \eqn{P(D|X=x)}. The options are:
\itemize{
\item \bold{``L'',} for \bold{L}inear logistic regression models.
\item \bold{``S'',} for \bold{S}mooth models.
}
}
}
\details{
\itemize{
\item If \code{meth} = \bold{``L''}, the logit transformation of the predicitive model is approximated by a linear logistic regression model:
\deqn{P (D|X=x) = 1/(1 + \exp{- \{ \beta_0 + \beta_1 x \}),}}
with \eqn{\beta_0, \beta_1 \in {\cal R}}.

\item If \code{meth} = \bold{``S''}, the logit transformation of the predicitive model is estimated by the smooth logistic regression,
\deqn{ P(D | X=x) = 1 / ( 1 + \exp \{ - s(x) \}),} being \eqn{s(\cdot)} the smooth function (splines, \doi{10.1002/sim.4780080504}).
}
}
\value{The returned value is a list with the two components:
\item{marker}{vector containing the ordered marker values.}
\item{probs}{vector with the probabilities corresponding to each marker value estimated through the predictive model.}
}



\seealso{
\code{sMS_binout} and \code{sMSROC}
}
