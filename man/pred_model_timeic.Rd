\name{pred_model_timeic-internal
}
\alias{pred_model_timeic
}
\title{Predictive model in prognosis scenarios (I)
}
\description{Estimation of the predictive model in prognosis scenarios under interval censorship.
}
\usage{pred_model_timeic(marker, left, right, outcome, time, meth)
}
\arguments{
\item{marker}{vector with the biomarker values.
}
\item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios and interval censorship and ignored in other situations.
}
\item{right}{vector with the upper edges of the observed intervals. It is mandatory as well in prognosis scenarios and interval censorship and ignored in other situations.
}
\item{outcome}{vector with the condition of the subjects as positive, negative or unknown at the considered time \code{time}.
}
\item{time}{point of time at which the sMS ROC curve estimator will be computed.
}
\item{meth}{ method for approximating the predictive model \eqn{P(D|X=x)}. The options are:
\itemize{
\item \bold{``L'',} for proportional hazards regression models taking into account the observation intervals (see Details).
\item \bold{``S'',} for proportional hazards regression models \bold{without} taking into account the observation intervals (see Details).
}
}
}
\details{
\itemize{
\item If \code{meth} = \bold{``L''}, the event times are assumed to come from a Cox proportional hazards regression model and the predictive model is estimated as indicated in
\doi{10.1080/00949655.2020.1736071}.
 \deqn{ P (T \leq t \;|\; X=x) = \frac{S(U|x) - S(t|x) }{S(U|x) - S(V|x)},} where \eqn{U = \min{\{t, L\}}} and \eqn{V = \max {\{t, R\}}}, being L and R the random variables that stand for the edges of the observable interval containing the event time.
\item If \code{meth} = \bold{``S''}, the approximation is done by \deqn{P (T\leq t \;|\; X=x) = 1 - S(t|x),} being \eqn{S(\cdot)} the survival function at time \bold{t} given the marker value, estimated through a proportional hazard model for interval censored data according to \doi{10.2307/2530698}.
}
}

\value{The returned value is a list with three components:

\item{marker}{vector containing the ordered marker values.
}
\item{probs}{vector with the probabilities corresponding to each marker value estimated through the predictive model.
}
\item{outcome}{vector with the condition of the subjects as positive, negative or censored at the considered time \code{time}.}

 }

\seealso{
\code{sMS_timeic} and \code{sMSROC}
}


