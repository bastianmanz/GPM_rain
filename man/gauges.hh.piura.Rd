\name{gauges.hh.piura}
\alias{gauges.hh.piura}
\docType{data}
\title{Example dataset of 30-min gauges in the Piura region, Peru.}
\description{Example dataset of 14 gauges in 30-min resolution in the Piura region of Peru. The dataset has a typical spread-sheet layout, such as obtained when reading a .csv, .xls or .txt file with coordinates, elevation and station names specified in the first rows and the temporal information thereafter.}
\usage{data("gauges.hh.piura")}
\format{
  A data frame with 4420 observations on the following 15 variables.
  \describe{
    \item{\code{V1}}{a factor with levels \code{Latitude} \code{Longitude} \code{Elevation} \code{DATE} in the first four rows and 30-min time-stamps from \code{2014-04-01 00:00:00} to \code{2015-04-01 00:00:00} thereafter.}
    \item{\code{V2}}{a factor with levels specifying the latitude, longitude, elevation and station name/ code in the first four rows (in that order) and thereafter containing the gauge measurements corresponding to the dates in V1. V3 to V15 are set up identifally to V2.}
    \item{\code{V3}}{a factor identical to V2.}
    \item{\code{V4}}{a factor identical to V2.}
    \item{\code{V5}}{a factor identical to V2.}
    \item{\code{V6}}{a factor identical to V2.}
    \item{\code{V7}}{a factor identical to V2.}
    \item{\code{V8}}{a factor identical to V2.}
    \item{\code{V9}}{a factor identical to V2.}
    \item{\code{V10}}{a factor identical to V2.}
    \item{\code{V11}}{a factor identical to V2.}
    \item{\code{V12}}{a factor identical to V2.}
    \item{\code{V13}}{a factor identical to V2.}
    \item{\code{V14}}{a factor identical to V2.}
    \item{\code{V15}}{a factor identical to V2.}
  }
}
\details{Example dataset of 14 gauges in 30-min resolution in the Piura region of Peru. The dataset has a typical spread-sheet layout, such as obtained when reading a .csv, .xls or .txt file with coordinates, elevation and station names specified in the first rows and the temporal information thereafter.}
\source{iMHEA monitoring programme (Wouter Buytaert, Imperial College London)}
\examples{
data(gauges.hh.piura)
gauges <- gauges.hh.piura
dim(gauges)
head(gauges)
}
\keyword{datasets}
