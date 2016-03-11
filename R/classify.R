#' Function to classify gauge time-series (zoo object) based on some prescribed aggregation level.
#'
#' @description This method classifies a gauge time-series based on a pre-defined classification ("aggregation_level"). Currently only a seasonal classification is implemented. The output is added to the gauge spatial object as an additional spatial data frame column.
#'
#' @param gauge.ts A zoo object containing the gauge observations.
#' @param aggregation_level A character string defining the classification criteria.
#'
#' @details If aggregation_level is "season", the date strings of the zoo object guage.ts are converted to seasonal indicators, i.e. 1 (DJF), 2 (MAM), 3 (JJA), 4 (SON).
#'
#' @return Object \code{ts.classification} A zoo object with a single data column defining the assigned class for each time-step.
#'
#' @examples
#' data(gauges) # STFDF object with gauge observations
#' gauge.sp <- gauges@sp
#' gauge.ts <- as(gauges[,,1],"xts")
#' colnames(gauge.ts) <- gauge.sp$estacion
#' ts.classification <- classify(gauge.ts,aggregation_level="season")
#' summary(ts.classification)

classify <- function(gauge.ts,aggregation_level="season"){

  if(aggregation_level=="season"){

    mtly_index <- as.numeric(format(index(gauge.ts),"%m"))

    season <- rep(NA,length(mtly_index))
    season[which(mtly_index<=2 | mtly_index==12)]=1 #DJF
    season[which(mtly_index>=3 & mtly_index<=5)]=2 # MAM
    season[which(mtly_index>=6 & mtly_index<=8)]=3 # JJA
    season[which(mtly_index>=9 & mtly_index<=11)]=4 # SON

    ts.classification <- as.xts(zoo(season,index(gauge.ts)))
    colnames(ts.classification)=aggregation_level

  } else { print("This classification has not yet been implemented!")}

  return(ts.classification)

}
