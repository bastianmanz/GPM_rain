% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bias_correct.R
\name{bias_correct}
\alias{bias_correct}
\title{Function to bias correct colocated gpm pixels against gauge observations}
\usage{
bias_correct(gauge.sp, gauge.ts, gpm.sp, gpm.ts)
}
\arguments{
\item{gauge.sp}{A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.}

\item{gauge.ts}{A zoo object containing the gauge observations.}

\item{gpm.sp}{A SpatialPointsDataFrame object containing the locations and numeric identifiers of the GPM-IMERG pixels.}

\item{gpm.ts}{A zoo object containing the gpm-imerg observations.}
}
\value{
Object \code{gpm.ts} A zoo object identical in dimension and format to the input gpm.ts with bias corrected values.
}
\description{
This method determines and corrects a multiplicative bias to the satellite data (gpm time-series).
}
\details{
At every time-step, the ratio of the mean of all active recording gauges to the mean of all corresponding, colocated gpm pixels within the domain is determined. This ratio is then applied to all (gauged and ungauged) GPM pixels at that time-step. Depending on the completeness of the gauge records the number of active recording gauges at a particular time-step may vary, which will impact in the accuracy of the results.
}
\examples{
data(gpm) # STFDF object with gpm imerg observations
data(gauges) # STFDF object with gauge observations
gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)
gauge.sp <- gauges@sp
gauge.ts <- as(gauges[,,1],"xts")
colnames(gauge.ts) <- gauge.sp$estacion
gpm.sp <- gpm@sp
gpm.ts <- as(gpm[,,1],"xts")
colnames(gpm.ts) <- gpm.sp$estacion
gauge.list <- combine_gauges(gauge.sp,gauge.ts)
gauge.sp <- gauge.list$sp
gauge.ts <- gauge.list$ts
gpm.bc.ts <- bias_correct(gauge.sp,gauge.ts,gpm.sp,gpm.ts)
}

