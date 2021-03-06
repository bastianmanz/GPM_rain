% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvgm.R
\name{cvgm}
\alias{cvgm}
\title{Function to fit climatological semi-variograms}
\usage{
cvgm(gauge.ts, gauge.sp, gpm.ts, gpm.sp, ts.classification, class_select = 1,
  rain_threshold = 1e-04, vgm.type = "KED", vgm.model = "Exp",
  class.dist = 5)
}
\arguments{
\item{gauge.ts}{A zoo object containing the gauge observations.}

\item{gauge.sp}{A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.}

\item{gpm.ts}{A zoo object containing the gpm-imerg observations.}

\item{gpm.sp}{A SpatialPointsDataFrame object containing the locations and numeric identifiers of the GPM-IMERG pixels.}

\item{ts.classification}{A zoo object with a single data column defining the assigned class for each time-step.}

\item{class_select}{A positive integer defining the column from the object \code{ts.classification} from which the climatological classification is to be obtained.}

\item{rain_threshold}{A numeric value defining the level below which all rainfall will be ignored and measurements treated as zeroes. This level can be set to the gauge detection limit, which typically is 0.1 or 0.2 mm/hr. However, the function requires a positive value to be supplied. Hence, if only actual zero-measurements are supposed to be treated as zeroes, use a miniscule value, such as 1e-4 (default).}

\item{vgm.type}{A character string defining which data the variogram is to be calculated from. This can be either ordinary kriging (OK), kriging with external drift (KED; default) or the satellite (SAT). In the case of OK, only the gauge data is used. In the case of KED, the satellite is used as a secondary variable when defining the semi-variogram. In the case of SAT the semi-variogram is calculated from the full grid of satellite observations. If subsequently using the semi-variograms for kriging, use the same vgm.type in cvgm and in DK.}

\item{vgm.model}{A character string defining the theoretical semi-variogram model that is to be fitted to the averaged and normalized experimental semi-variogram. These models are derived from the gstat package. Type vgm() to obtain a full list. Multiple models can be provided and the one with the lowest sum of squared errors is selected and returned.}

\item{class.dist}{A numeric value defining the width (horizontal distance) in kilometres over which the experimental variograms within a climatological class are binned. Increasing this value will smooth out the variability and improve the fitting of the theoretical variogram, but may lead to numerical artefacts. A low value typically results in a lot of noise. The value should not be larger than the grid cell size fo the satellite product/ target grid (0.1 degree in the case of GPM-IMERG).}
}
\value{
Object \code{clim_vgm} A list containing four entries for each, occurrence and amounts. An event identifier is provided that defines the climatological class for event time-step/ event. The second entry contains the fitted models required for the kriging. The third entry is a data frame summarizing the climatological semi-variogram and number of events for each climatological class. The last entry contains the semi-variogram plots. Type names(clim_vgm) to obtain the names of all entries.
}
\description{
This function undergoes a number of steps to calculate experimental (empirical) semi-variograms of gauge-based time-series, aggregates these based on the previously defined climatological classification approach, and, ultimately, following a normalization step, fits theoretical semi-variogram models for each class. The motivation for claculating climatological variograms lies in the fact that, depending on data availability, fitting theoretical variogram models at each individual time-step may not be successful, making climatological variograms a more robust approach. This step is essential to for Double Kriging.
}
\details{
The climatological variogram function is a required pre-requisite for the Double Kriging (DK) step as it defines the spatial variance structure (semi-variogram) that is required for prediction in DK. The input time-series data is converted into an occurrence time-series (binary indicator) or a normal-score transformed Gaussian rainfall intensity time-series. For both of these time-series an experimental semi-variogram is determined at each time-step. After this is done for the full time-series, the semi-variograms are classified depending on the pre-defined climatological classification and normalized and averaged for each class.
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
ts.classification <- classify(gauge.ts,aggregation_level="season")
clim.vgms.ked <- cvgm(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification,class_select=1,rain_threshold=0.001,vgm.type="KED",vgm.model=c("Exp","Gau","Sph","Mat"), class.dist=5)
class(clim.vgms.ked)
names(clim.vgms.ked)
clim.vgms.ked$occ_vgm_results
clim.vgms.ked$amts_vgm_results
print(clim.vgms.ked$occ_plot[[1]], split=c(1,1,4,2), more=TRUE)
print(clim.vgms.ked$occ_plot[[2]], split=c(2,1,4,2), more=TRUE)
print(clim.vgms.ked$occ_plot[[3]], split=c(3,1,4,2), more=TRUE)
print(clim.vgms.ked$occ_plot[[4]], split=c(4,1,4,2), more=TRUE)
print(clim.vgms.ked$amts_plot[[1]], split=c(1,2,4,2), more=TRUE)
print(clim.vgms.ked$amts_plot[[2]], split=c(2,2,4,2), more=TRUE)
print(clim.vgms.ked$amts_plot[[3]], split=c(3,2,4,2), more=TRUE)
print(clim.vgms.ked$amts_plot[[4]], split=c(4,2,4,2), more=TRUE)
}

