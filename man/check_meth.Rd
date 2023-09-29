\name{check_meth-internal}
\alias{check_meth}
\title{Check the method for estimating the predictive model}
\description{When the predictive model is entered manually by the user, this function ensures that no other method by default is used to compute it.}
\usage{check_meth(meth, probs)}
\arguments{
  \item{meth}{ method for approximating the predictive model \eqn{P (D | X=x)}.}
  \item{probs}{ vector containing the probabilities corresponding to the predictive model when it is entered manually by the user.}
}
\details{
If the predicitve model has been manually indicated, this function sets the parameter to {\strong{"M"}} ignoring other options. In this way, none of the function computing the predictive model will be called.
}
\value{A list with one component:
\item{meth}{ value set up for the method parameter.It can take eiter the value entered by the user or its default, if the predictive model was not manually indicated, or the value {\strong{"M"}}, when the predicitive model was entered in the parameter {\code{probs}}.
}
}

