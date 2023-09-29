\name{diabet}
\alias{diabet}
\docType{data}
\title{Diabetes dataset
}
\description{This dataset contains part of the Diabetes Dataset (see References), courtesy of Dr John Schorling from the Department of Medicine, University of Virginia School of Medicine.
This version contains 3 variables on 403 subjects interviewed to understand the prevalence of several cardiovascular risks factors in central Virginia for African Americans.
}
\usage{data("diabet")}
\format{
  A data frame with 403 observations on the following 3 variables.
  \describe{
    \item{\code{stab.glu}}{a numeric vector indicating the level of stabilized glucose.}
    \item{\code{glyhb}}{a numeric vector indicating the level of glycosolated hemoglobin.}
    \item{\code{age}}{age in years of the participants.}
    \item{\code{diab}}{a numeric vector indicating whether the subject is diagnosed as diabetic (value = 1) or not (value = 0). }
  }
}
\details{The \bold{diab} variable is not present in the original dataset. Here, values of glycosolated hemoglobin > 7.0 were taken as a positive diagnosis of diabetes (\bold{diab} = 1) and those of glycosolated hemoglobin <= 7.0 as a negative diagnosis (\bold{diab} = 0).
}
\source{Full dataset can be downloaded at \url{https://hbiostat.org/data}.
}
\references{Willems JP, Saunders JT, Hunt DE, Schorling JB. Prevalence of coronary heart disease risk factors among rural blacks: a community-based study. South Med J. 1997 Aug;90(8):814-20. PMID: 9258308.}
\examples{
data(diabet)
summary(diabet)
}

