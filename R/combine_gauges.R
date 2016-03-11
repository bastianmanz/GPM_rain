#' Function to combine gauges within grid cells for kriging
#'
#' @description This function identifies if multiple gauges are located in the same grid cell (e.g. GPM pixel) and averages their time-series to obtain a single best ground-based estimate to compare and merge with the satellite data.
#'
#' @param gauge.ts A zoo object containing the gauge observations.
#' @param gauge.sp A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.
#'
#' @details This function averages (simple mean) multiple gauges if these are located in the same grid cell (i.e. GPM pixel) such that a single time-series for the ground-based estimate is obtained for every grid node. In order for this function to work, the function \code{colocate} needs to be run previously such that the required GPM pixel identifiers are assigned to each gauge location.
#'
#' @return Object \code{gauge.list} A list with two entries: the gauge spatial data identical to that in the object \code{gauge.sp} and a revised time-series object with a single unique time-series corresponding to each GPM grid node that contained gauges.
#'
#' @examples
#' data(gpm) # STFDF object with gpm imerg observations
#' data(gauges) # STFDF object with gauge observations
#' gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)
#' gauge.sp <- gauges@sp
#' gauge.ts <- as(gauges[,,1],"xts")
#' colnames(gauge.ts) <- gauge.sp$estacion
#' gauge.list <- combine_gauges(gauge.sp,gauge.ts)
#' gauge.sp <- gauge.list$sp
#' gauge.ts <- gauge.list$ts

combine_gauges <- function(gauge.sp,gauge.ts){

  duplicates <- gauge.sp@data$GPM_PIXEL[which(duplicated(gauge.sp@data$GPM_PIXEL))]

  for(i in 1:length(duplicates)){

    id_i <- which(gauge.sp@data$GPM_PIXEL==duplicates[i])

    gauge.sp <- gauge.sp[-id_i[2:length(id_i)],]

    gauge.ts[,id_i[1]] <- rowMeans(gauge.ts[,id_i],na.rm=T)
    gauge.ts <- gauge.ts[,-id_i[2:length(id_i)]]

  }

  gauge.list <- list()
  gauge.list$sp <- gauge.sp
  gauge.list$ts <- gauge.ts
  return(gauge.list)
}
