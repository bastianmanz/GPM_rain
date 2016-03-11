\name{gpm}
\alias{gpm}
\docType{data}
\title{Example Daily-Aggregated GPM-IMERG Data}
\description{Example GPM-IMERG Production version 3 data from the Lima region, Peru. The 30-min observations are aggregated to the daily level. The data is stored in a STFDF object (spacetime) that contains spatial information on the GPM pixel grid nodes and temporal information on the measurement time-series at each grid node.}
\usage{data("gpm")}
\format{A STFDF (spacetime) object with spatial information on the GPM grid node locations and temporal information (i.e.time-series of the GPM pixels.)}
\details{
Example GPM-IMERG Production version 3 data from the Lima region, Peru. The 30-min observations are aggregated to the daily level. }
\source{NASA PPS http://pmm.nasa.gov/data-access/downloads/gpm}
\examples{
data(gpm) # STFDF object with gpm imerg observations
str(gpm)
gpm.sp <- gpm@sp
gpm.ts <- as(gpm[,,1],"xts")
colnames(gpm.ts) <- gpm.sp$estacion
summary(gpm.sp)
class(gpm.ts)
dim(gpm.ts)
}
\keyword{datasets}
