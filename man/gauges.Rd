\name{gauges}
\alias{gauges}
\docType{data}
\title{Example dataset of daily gauges in the region of Lima, Peru.}
\description{Example dataset of daily gauges in the region of Lima, Peru. The dataset is saved in STFDF (spacetime) format and allows for extracting spatial information and daily time-series (zoo).}
\usage{data("gauges")}
\format{A STFDF (spacetime) object with spatial information on the gauge locations and temporal information (i.e.time-series of the gauges.)}
\details{Example dataset of daily gauges in the region of Lima, Peru.}
\source{SENAMHI, Lima, Peru}
\examples{
data(gauges)
str(gauges)
data(gpm) # STFDF object with gpm imerg observations
gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)
gauge.sp <- gauges@sp
gauge.ts <- as(gauges[,,1],"xts")
colnames(gauge.ts) <- gauge.sp$estacion
summary(gauge.sp)
summary(gauge.ts)
}
\keyword{datasets}
