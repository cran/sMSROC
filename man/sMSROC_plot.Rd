\name{sMSROC_plot}
\alias{sMSROC_plot}
\title{Plot of the sMS ROC curve estimate
}
\description{Provides informative plots of the sMS ROC curve estimates.
}
\usage{
sMSROC_plot(sMS, m.value)
}
\arguments{
  \item{sMS}{object of class \code{sMS} returned from function \code{\link{sMSROC}}.
}
  \item{m.value}{marker value. It is an optional parameter that, when indicated, adds over the graphic of the ROC curve, the point which corresponds to that marker value.
}
}
\details{
The function provides two types of graphics:
\itemize{
\item{A basic plot approximating the ROC curve by the pairs given by the sequences {\bold{1 - SP}} and {\bold{SE}}, from the \code{sMSROC} object. The layers \code{geom_roc()} and \code{roc_style()} from the \code{\link{plotROC}} package  were added to this plot, which make possible to take advantage of the functionality of this package.
}

\item{A customized graphic of the ROC curve whose class is \code{ggplot}, obtained approximating the sequences {\bold{1 - SP}} and {\bold{SE}}. When te parameter \code{m.value} is indicated, the final plot displays over the ROC curve estimate the point that corresponds to the entered value.
}
}
}
\value{A list with the following elements:
\item{basic.plot}{object that can be used and customized by the tools from the \code{\link{plotROC}}  package.
}
\item{roc.plot}{object of class  \code{ggplot}. Although it is already customized (title, colors, axis labels, ..., etc.) the end-users can make their own changes by adding the corresponding layers, with the available tools from the \code{ggplot2} package.
}
}
\seealso{
\code{sMSROC}
}
\examples{
#  Example of the use of the plot.sMSROC function
data(ktfs)
DT = ktfs
ROC <- sMSROC(marker = DT$score,
              status = DT$failure,
              observed.time = DT$time,
              time = 5,
              meth = "S")
plot <- sMSROC_plot(sMS = ROC, m.value = 4.2)
plot$basicplot; plot$rocplot
}
