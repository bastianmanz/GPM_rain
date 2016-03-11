#' Function to bias correct colocated gpm pixels against gauge observations
#'
#' @description This method determines and corrects a multiplicative bias to the satellite data (gpm time-series).
#'
#' @param gauge.ts A zoo object containing the gauge observations.
#' @param gauge.sp A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.
#' @param gpm.ts A zoo object containing the gpm-imerg observations.
#' @param gpm.sp A SpatialPointsDataFrame object containing the locations and numeric identifiers of the GPM-IMERG pixels.
#'
#' @details At every time-step, the ratio of the mean of all active recording gauges to the mean of all corresponding, colocated gpm pixels within the domain is determined. This ratio is then applied to all (gauged and ungauged) GPM pixels at that time-step. Depending on the completeness of the gauge records the number of active recording gauges at a particular time-step may vary, which will impact in the accuracy of the results.
#'
#' @return Object \code{gpm.ts} A zoo object identical in dimension and format to the input gpm.ts with bias corrected values.
#'
#' @examples
#' data(gpm) # STFDF object with gpm imerg observations
#' data(gauges) # STFDF object with gauge observations
#' gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)
#' gauge.sp <- gauges@sp
#' gauge.ts <- as(gauges[,,1],"xts")
#' colnames(gauge.ts) <- gauge.sp$estacion
#' gpm.sp <- gpm@sp
#' gpm.ts <- as(gpm[,,1],"xts")
#' colnames(gpm.ts) <- gpm.sp$estacion
#' gauge.list <- combine_gauges(gauge.sp,gauge.ts)
#' gauge.sp <- gauge.list$sp
#' gauge.ts <- gauge.list$ts
#' gpm.bc.ts <- bias_correct(gauge.sp,gauge.ts,gpm.sp,gpm.ts)

bias_correct <- function(gauge.sp,gauge.ts,gpm.sp,gpm.ts){

  # calculate ts of multiplative bias ratio
  mn.gpm.coloc.ts <- rowMeans(gpm.ts[,unique(gauge.sp@data$GPM_PIXEL)],na.rm=T)
  mn.gauge.ts <- rowMeans(gauge.ts,na.rm=T)

  ratio.ts <- mn.gauge.ts/mn.gpm.coloc.ts
  ratio.ts[which(!is.finite(ratio.ts))]=NA # eliminate ts where both = 0

  # apply bias ratio
  gpm.ts[which(!is.na(ratio.ts)),] <- gpm.ts[which(!is.na(ratio.ts)),] * ratio.ts[which(!is.na(ratio.ts))]

  return(gpm.ts)

}
