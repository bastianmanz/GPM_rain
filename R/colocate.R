#' Function to colocate gauges with corresponding GPM pixel
#'
#' @description This function finds the closest ("colocated") GPM-IMERG pixel within a pre-defined radius to the provided gauge locations. This allows for comparing corresponding gauge and satellite time-series.
#'
#' @param gauges A STFDF object containing the gauge spatial information and time-series data.
#' @param gpm a STFDF object containing the GPM-IMERG spatial information and time-series data.
#' @param resolution An integer describing the satellite spatial resolution.
#' @param longlat A logical object indicating if object is in long/lat or not. TRUE by default.
#'
#' @details The closest GPM-IMERG pixel to each gauge is found by iteratively calculating the Euclidian (if longlat=FALSE) or Great Circle distance between each gauge location and all GPM-IMERG pixel locations using the spDistsN1() function from the sp package. The closest GPM pixel is then identified as the colocated pixel. The parameter "resolution" defines the satellite grid resolution, which is used to define the allowable tolerance limit for the distance between the gauge and the colocated GPM pixel use the pythagoras theorem.
#'
#' @return Object \code{gauges} A STFDF object identical to the input object with an additional spatial data column ("GPM_PIXEL" that defines the ID of the colocated GPM pixel for each gauge location.
#'
#' @examples
#' data(gpm) # STFDF object with gpm imerg observations
#' data(gauges) # STFDF object with gauge observations
#' gauges <- colocate(gauges,gpm,resolution=0.1,longlat=TRUE)

colocate <- function(gauges, gpm, resolution,longlat=TRUE){

  pts <- gauges@sp
  coords <- gpm@sp

  # get location of gauges in gpm pixels
  loc1 <- numeric()
  loc_dists <- numeric()

  for (i in 1:length(pts)) {
    loc1[i] <- which.min(spDistsN1(coords,pts[i,],longlat))
    loc_dists[i] <- min(spDistsN1(coords,pts[i,] ,longlat))
  }

  # remove gauges which have no pixel within length of pixel radius
  if(longlat==TRUE){
    x=sqrt(2*((resolution*100)^2)) # pythagoras, convert from deg latlong to km
  } else{
    x=sqrt(2*(resolution^2)) # pythagoras
  }
  loc <- loc1[c(which(loc_dists < x))]

  gauges@sp@data <- cbind(gauges@sp@data,loc)
  colnames(gauges@sp@data)[ncol(gauges@sp)] <- "GPM_PIXEL"

  return(gauges)
}
