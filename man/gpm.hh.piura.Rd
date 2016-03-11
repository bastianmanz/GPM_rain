\name{gpm.hh.piura}
\alias{gpm.hh.piura}
\docType{data}
\title{Example 30-min GPM-IMERG Data}
\description{Example GPM-IMERG Production version 3 data from the Piura region, Peru at the original 30-min resolution (UTC time-zone).The data is stored in a STFDF object (spacetime) that contains spatial information on the GPM pixel grid nodes over a 1.0 x 1.0 degree space (100 nodes) and temporal information on the measurement time-series at each grid node.}
\usage{data("gpm.hh.piura")}
\format{A STFDF (spacetime) object with spatial information on the GPM grid node locations and temporal information (i.e.time-series of the GPM pixels.)}
\details{Example GPM-IMERG Production version 3 data from the Piura region, Peru at the original 30-min resolution (UTC time-zone) and 0.1 degree spatial resolution.}
\source{NASA PPS http://pmm.nasa.gov/data-access/downloads/gpm}
\examples{
data(gpm.hh.piura) # STFDF object with gpm imerg observations
gpm <- gpm.hh.piura
str(gpm)
gpm.sp <- gpm@sp
gpm.ts <- as(gpm[,,1],"xts")
colnames(gpm.ts) <- gpm.sp$estacion
summary(gpm.sp)
class(gpm.ts)
dim(gpm.ts)}
\keyword{datasets}
