\name{fibrosis}
\alias{fibrosis}
\docType{data}
\title{Fibrosis dataset
}
\description{Synthetic dataset generated to simulate data from a study that aimed to assess the predictive ability of a constructed score to determine the worsening in the fibrosis stage in individuals infected by the hepatitis C (HC) virus. When participants underwent a target revision, their fibrosis stage, certain polymorphisms, and other clinical variables were collected. Highest stages of fibrosis were considered a worsening in the disease. See References for more information about the study.

}
\usage{data("fibrosis")}
\format{
  A data frame with 722 observations and the following variables:
  \describe{
    \item{\code{Id}}{Identification label for each participant.
    }
    \item{\code{Score}}{Score proposed to stratify participants infected by the HC virus according to the risk of a worsening in their fibrosis stage.
    }
    \item{\code{Start}}{Lower edge of the observable interval.
    }
    \item{\code{Stop}}{Upper edge of the observable interval. It can take the value infinity, represented as Inf.
    }
  }
}
\source{Synthetic dataset.
}
\references{Vidal-Castineira JR. al. Genetic contribution of endoplasmic reticulum aminopeptidase 1 polymorphisms to liver fibrosis progression in patients with HCV infection. Journal of Molecular Medicine, 98:1245-1254, 2020. \doi{10.1007/s00109-020-01948-1}
}
\examples{
data(fibrosis)
summary(fibrosis)
}

