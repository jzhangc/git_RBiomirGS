% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{rbiomirgs_bar}
\alias{rbiomirgs_bar}
\title{rbiomirgs_bar}
\usage{
rbiomirgs_bar(
  gsadfm,
  export_name = NULL,
  gs.name = FALSE,
  n = "all",
  signif_only = FALSE,
  p_adj = TRUE,
  q_value = 0.05,
  y.rightside = FALSE,
  Title = NULL,
  xLabel = "Log odds ratio",
  yLabel = NULL,
  errorbarWidth = 0.2,
  xTxtSize = 10,
  yTxtSize = 10,
  plotWidth = 250,
  plotHeight = 230
)
}
\arguments{
\item{gsadfm}{Input dataframe. The dataframe is the output from \code{\link{rbiomirgs_logistic}}.}

\item{export_name}{String. File name prefix for output plots. Default is \code{NULL}.}

\item{gs.name}{Whether or not to display gene set names on the y-axis. Default is \code{FALSE}.}

\item{n}{Number of gene sets to plot. Options are \code{"all"} and an integer number. Default is \code{"all"}.}

\item{signif_only}{Whether to plot only the significantly enriched GS. Default is \code{FALSE}.}

\item{p_adj}{Set when \code{signif_only = TRUE}, whether to use adjusted p value for thresholding. Default is \code{TRUE}.}

\item{q_value}{Set when \code{signif_only = TRUE}, the significance threshold. Default is \code{0.05}.}

\item{y.rightside}{Whether or not to display right side y-axis. Default is \code{FALSE}.}

\item{Title}{Title of the plot. Default is \code{NULL}.}

\item{xLabel}{X-axis label. Default is \code{"Log  odds ratio"}.}

\item{yLabel}{Y-axis label. Default is \code{NULL}.}

\item{plotWidth}{Set the width of the plot. Default is \code{250}.}

\item{plotHeight}{Set the height of the plot. Default is \code{230}.}
}
\value{
Outputs a \code{pdf} boxplot figure.
}
\description{
Bar graph function
}
\details{
The function produces either a bar graph or a vertical version. When \code{gs.name = FALSE}, the function produces a plot with x-axis as gene set. When \code{gs.name = TRUE}, the plot will have gene set names displayed on the y-axis, and x-axis will be used for log oods ratio. The error bar is standard error of the coefficient (log odds ratio).
}
\examples{
\dontrun{

# horizontal version
rbiomirgs_bar(gsadfm = mirna_mrna_GS, n = "all", y.rideside = TRUE, yTxtSize = 8, plotWidth = 300, plotHeight = 200, xLabel = "Gene set", yLabel = "Log odds ratio")

# vertical version
 rbiomirgs_bar(gsadfm = mirna_mrna_GS, n = 60, yAxis = TRUE,
 yTxtSize = 8, plotWidth = 300, plotHeight = 200, xLabel = "Log odds ratio")

}
}
