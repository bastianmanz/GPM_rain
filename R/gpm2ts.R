#' Function to extract xts from STFDF for GPM-IMERG data
#'
#' @description This function is a computationally efficient way to convert a large-scale GPM-IMERG STFDF (spacetime) object into an xts object, while also truncating the temporal extent to that of the supplied gauge time-series (xts object).
#'
#' @param gpm A STFDF (spacetime) object containing the GPM-IMERG spatial information (i.e. grid nodes with attributes) and temporal information (i.e. precipitation time-series for each grid node location).
#' @param gauge.ts A zoo object containing the gauge observations.
#'
#' @details This function is a computationally efficient way to convert a large-scale GPM-IMERG STFDF (spacetime) object into an xts object, while also truncating the temporal extent to that of the supplied gauge time-series (xts object).
#'
#' @return Object \code{gpm.ts} A zoo, xts object with temporal extent from the start date to the end date of the corresponding gauge zoo object and GPM-IMERG estimates in columns for the GPM grid nodes. To ensure that the temporal extent of gpm.ts equals that of gauge.ts, ensure that gauge.ts has the same time-step length as GPM-IMERG (30-min) and time-stamps coincide.
#'
#' @examples
#' data(gpm.hh.piura) # STFDF object with gpm imerg observations for the Piura region in Peru
#' gpm <- gpm.hh.piura
#' data(gauges.hh.piura) # STFDF object with gauge observations for the Piura region in Peru
#' gauges <- gauges.hh.piura
#' gauges <- spacetime_gauge(gauges,projection="+proj=longlat +ellps=WGS84 +datum=WGS84", timezone=-5)
#' gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)
#' gauge.sp <- gauges@sp
#' gpm.sp <- gpm@sp
#' gauge.ts <- as(gauges[,,1],"xts")
#' colnames(gauge.ts) <- gauge.sp$estacion
#' gpm.ts <- gpm2ts(gpm,gauge.ts)

gpm2ts <- function(gpm,gauge.ts){

  gpm_data <- matrix(data=gpm@data[,1],ncol=length(gpm@sp),nrow=length(gpm@time),byrow=TRUE)
  gpm.ts <- as.xts(zoo(gpm_data,index(gpm@time)))

  start <- which(index(gpm.ts)==min(index(gauge.ts)))
  end <- which(index(gpm.ts)==max(index(gauge.ts)))

  gpm.ts <- gpm.ts[start:end,]

  return(gpm.ts)

}
