\name{sMSROC}
\alias{sMSROC}
\title{sMS ROC curve estimator computation}
\description{Core function for computing the sMS ROC estimator which fits the estimation of the ROC curve when the outcome of interest is time-dependent (\bold{prognosis} scenarios) and when it is not (\bold{diagnosis} scenarios).}
\usage{
sMSROC(marker, status, observed.time, left, right, time,
       meth, grid, probs, sd.probs,
       conf.int, ci.cl, ci.meth, ci.nboots, parallel, ncpus, all)
}
\arguments{
  \item{marker}{vector with the biomarker values.}

  \item{status}{numeric response vector. The highest value is assumed to stand for the subjects having the event under study. The lowest one, for those who do not. Any other value will not be considered. It is a mandatory parameter in diagnosis scenarios.}

  \item{observed.time}{vector with the observed times for each subject, for prognosis scenarios under right censorship. Notice that these values may be the event times or the censoring times.}

 \item{left}{vector containing the lower edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations.}

  \item{right}{vector with the upper edges of the observed intervals. It is mandatory in prognosis scenarios under interval censorship and ignored in other situations. The infinity is admissible as value (indicated as \bold{inf}).}

  \item{time}{ point of time at which the sMS ROC curve estimator will be computed. The default value is 1.}

  \item{meth}{ method for approximating the predictive model \eqn{P(D|X=x)}. There are several options available:

\itemize{
\item \bold{``E'',} allocates to each individual their own condition as positive or negative. Those whose condition is unknown at time \code{time} are dismissed.

\item \bold{``L'',} for \bold{L}inear logistic regression and proportional hazards regression models (see Details).

\item \bold{``S'',} for \bold{S}mooth models (see Details).
}
}

  \item{probs}{ vector containing the probabilities corresponding to the predictive model when it has been externally computed. Only values within [0,1] are admissible.}

  \item{sd.probs}{ vector with the standard deviations of the probabilities entered in \code{probs}. It is an optional parameter.}

  \item{grid}{grid size for computing the AUC. Default value 1000.}
  \item{conf.int}{indicates whethet a conficence interval for the AUC will be computed (\bold{``T''}) or not (\bold{``F''}). The default value is (\bold{``F''}).}

  \item{ci.cl}{confidence level at which the confidence interval for the AUC will be provided. The default value is 95\%. This parameter is ignored when \code{conf.int} is set to \bold{``F''}.
  }
\item{ci.meth}{method for computing the confidence interval for the AUC. There are three options:
\itemize{
   \item \bold{``E''}, for the \bold{E}mpirical variance estimation.
   \item \bold{``V''}, for the theoretical \bold{V}ariance estimation.
   \item \bold{``B''}, for the \bold{B}ootstrap percentile approximation.

  The empirical method \bold{E} is taken as default value and the parameter is ignored too when \code{ci.cl} value is \bold{``F''}.
}
}
\item{ci.nboots}{number of boostrap samples to be run when \bold{B}oostrap is set as \code{ci.meth} parameter. The default value is 500 and it is not taken into account when no confidence interval is computed.}
\item{parallel}{indicates whether parallel computing will be done (\bold{``T''}) or not (\bold{``F''}) when computing the variance of the AUC through the methods \bold{``V''} and \bold{``B''}.
}
\item{ncpus}{number of CPUS that will be used when parallel computing is chosen. The default value is 1 and the maximum is 2.
}
\item{all}{parameter indicating whether all probabilities given by the predictive model should be considered (value \bold{``T''}) or just those corresponding to individuals whose condition as positive or negative is unknown (\bold{``F''}). The default value is (\bold{``T''}).
}
}

\details{The Two-stages mixed-subjects (sMSROC) ROC curve estimator links diagnosis and prognosis scenarios through a general predictive model (first
stage) and the weighted empirical estimator of the cumulative distribution function of the biomarker
(second stage).

The predictive model \eqn{P(D|X=x)} depicts the relationship between the biomarker and the binary response variable. It is approximated through the most suitable probabilistic model.

For \bold{diagnosis} scenarios:
\itemize{
\item If \code{meth} = \bold{``L''}, the logit transformation of the predicitive model is approximated by a linear logistic regression model:
\deqn{P (D|X=x) = 1/(1 + \exp{- \{ \beta_0 + \beta_1 x \}),}}
with \eqn{\beta_0, \beta_1 \in {\cal R}}.
\item If \code{meth} = \bold{``S''}, the logit transformation of the predicitive model is estimated by the smooth logistic regression,
\deqn{ P(D | X=x) = 1 / ( 1 + \exp \{ - s(x) \}),} being \eqn{s(\cdot)} the smooth function (splines, \doi{10.1002/sim.4780080504}).
}

Notice that the predictive model allows to compute the probability of being positive/negative even when the actual belonging group is unknown.

For \bold{prognosis} scenarios and right censorship:
\itemize{
\item If \code{meth} = \bold{``L''}, the event times are assumed to come from a Cox proportional hazards regression model:
 \deqn{ P (T \leq t \;|\; X=x) = 1 - \exp \{ - \Delta_0(t) \cdot \exp \{ \beta_0 + \beta_1 \cdot \log(x)\}\},} where \eqn{\Delta_0(\cdot)} is the baseline hazard function and \eqn{\beta_0, \beta_1 \in {\cal R}}.

\item If \code{meth} = \bold{``S''}, the approximation is done by \deqn{P (T\leq t \;|\; X=x) = 1 - \exp \{ - \Delta_0(t) \cdot \exp \{ s(x)\}\} } being \eqn{s(\cdot)} the smooth function (penalized splines, \doi{10.1111/1467-9868.00125}).}

Finally, for \bold{prognosis} scenarios and interval censorship:
\itemize{
\item{If \code{meth} = \bold{``L''}, the event times are assumed to come from a Cox proportional hazards regression model and the predictive model is estimated as indicated in \doi{10.1080/00949655.2020.1736071}.
 \deqn{ P (T \leq t \;|\; X=x) = \frac{S(U|x) - S(t|x) }{S(U|x) - S(V|x)},} where \eqn{U = \min{\{t, L\}}} and \eqn{V = \max {\{t, R\}}}, being L and R the random variables that stand for the edges of the observable interval containing the event time.
}

\item{If \code{meth} = \bold{``S''}, the approximation is done by \deqn{P (T\leq t \;|\; X=x) = 1 - S(t|x),} being \eqn{S(\cdot)} the survival function at time \bold{t} given the marker value, estimated through a proportional hazard model for interval censored data according to \doi{10.2307/2530698}.
}
}

The confidence intervals for the AUC can be computed in three different ways according to parameter \code{ci.meth}. When it is set to \bold{"E"} the variance of the AUC is estimated by the empirical
procedure and when the chosen option is \bold{"V"}, the theoretical approximation is used (see \doi{10.1515/ijb-2019-0097}). The third option in by using the Bootstrap percentile.

}

\value{The ouput is an objetc of class \code{sMSROC} with the following components:
\item{thres}{vector containing the biomarker values for which sensitivity and specificity were computed.}
\item{SE}{vector with the estimates of the sensitivity.}
\item{SP}{vector with the estimates of the specificity.}
\item{probs}{vector with the probabilities corresponding to the predictive model.}
\item{u}{vector containing the points between 0 and 1 at which the ROC curve estimator will be computed. Its size is determined by the \code{grid} parameter.}
\item{ROC}{ROC curve approximated at each point of the vector \code{u}.}
\item{auc}{ area under sMSROC curve estimator.}
\item{auc.ci.l}{lower edge of the confidence interval for the AUC.}
\item{auc.ci.u}{upper edge of the confidence interval for the AUC.}
\item{ci.cl}{confidence level at which the confidence interval for the AUC were computed.}
\item{ci.meth}{method chosen for computing the confidence interval for the AUC.}
\item{time}{point of time at which the sMS ROC curve estimator was computed in prognosis scenarios.}
\item{data}{list contaning several parameters used in the internal functions, when applicable:
\itemize{
\item data_type - type of scenario handled (diagnosis/prognosis, under right or interval censorship).
\item grid - grid size.
\item marker - vector with the biomarker values.
\item outcome  - vector with the condition of the individuals at time \code{time} as positive, negative or unknown.
\item ncpus - CPUs used if parallel computing was performed.
\item ci.nboots - number of bootstrap samples generated for computing the confidence intervals for the AUC.
\item parallel - was parallel computing performed?
\item meth - method used to compute the predictive model.
\item status - response vector.
\item observed.time - vector with the observed times for each subject.
\item left - vector with the lower edges of the observed intervals.
\item right - vector with the upper edges of the observed intervals.
}
}
\item{message}{table containing the warning messages generated during the execution of the function.}
}

\examples{
data(ktfs)
DT <- ktfs
sROC <- sMSROC(marker = DT$score, status = DT$failure,
               observed.time = DT$time, time = 5, meth = "L", conf.int = "T",
               ci.cl =0.90, ci.meth = "E")
}

\references{
{S. Díaz-Coto, P. Martínez-Camblor, and N. O. Corral-Blanco. Cumulative/dynamic ROC curve
estimation under interval censorship. Journal of Statistical Computation and Simulation, 90(9):1570–
1590, 2020. \doi{10.1080/00949655.2020.1736071}.
}


{
S. Díaz-Coto, N. O. Corral-Blanco, and P. Martínez-Camblor. Two-stage receiver operating-characteristic
curve estimator for cohort studies. The International Journal of Biostatistics, 17:117–137, 2021. \doi{10.1515/ijb-2019-0097}.
}

{Finkelstein, Dianne M. A Proportional Hazards Model for Interval-Censored Failure Time Data. Biometrics 42, no. 4 (1986): 845–54. \doi{10.2307/2530698}.
}

{Durrleman S, Simon R. Flexible regression models with cubic splines. Statistics in Medicine 1989; 8(5): 551-561.\doi{10.1002/sim.4780080504}
}

{Hurvich C, Simonoff J, Tsai CL. Smoothing parameter selection in nonparametric regression using an improved Akaike 1998. J.R. Statist. Soc. 60 271-293. \doi{10.1111/1467-9868.00125}
}

{B. Efron and R. J. Tibshirani. An Introduction to the Bootstrap. CRC press, 1994.}

}

