\name{ktfs}
\alias{ktfs}
\docType{data}
\title{KTFS dataset
}
\description{Dataset originally delivered in the \code{RISCA} package. It contains data from kidney transplant recipients for whom the Kidney Transplant Failure Score (KTFS) was collected. The KTFS is a score proposed by Foucher et al. (2010) (see References) to assess the recipients according to their risk of returning in dialysis.
}
\usage{data("ktfs")}
\format{
  A data frame with 2169 observations and the following 3 variables:
  \describe{
    \item{\code{time}}{a numeric vector depicting the follow-up time in years.
    }
    \item{\code{failure}}{a numeric vector indicating the graft failure at the end of the follow-up (1-Yes, 0-Censoring).
    }
    \item{\code{score}}{a numeric vector representing the  KTFS value.
    }
  }
}

\source{This dataset is available at \code{RISCA} package. More information about the KTFS score can be found at \url{https://www.divat.fr}.
}
\references{Foucher Y. al. A clinical scoring system highly predictive of long-term kidney graft survival. Kidney International, 78:1288-94, 2020. \doi{10.1038/ki.2010.232}.
}
\examples{
data(ktfs)
summary(ktfs)
}

