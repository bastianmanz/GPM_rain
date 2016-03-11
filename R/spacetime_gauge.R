#' Function to convert gauge information from a data.frame to STFDF (spacetime)
#'
#' @description This function converts gauge information including names, codes, coordinates and time-series from a typical spread-sheet layout (in .csv, .xls, or .txt) to a STFDF (spacetime) object. The time-stamps are also converted from the local time-zone to UTC. Use data(gauges.hh.piura) to obtain an example of the format of the data-frame for which this function is designed.
#'
#' @param gauges A dataframe containing gauge information in a typical format. Type data (gauges.piura) to see an example.
#' @param projection A character string detailing the projection CRS. See the \code{sp} package for details on this, in particular the function \code{proj4string} therein. The default(\code{projection}="+proj=longlat +ellps=WGS84 +datum=WGS84") is suitable if coordinates are in latitude/ longitude, using a WGS84 datum.
#' @param timezone a numeric value, describing the timezone difference to UTC (all GPM-IMERG data is in UTC) so that the gauge time-series can be converted to UTC. Negative value imply western hemisphere, e.g. Peru = -5.
#'
#' @details This function converts gauge information including names, codes, coordinates and time-series from a typical spread-sheet layout (in .csv, .xls, or .txt) to a STFDF (spacetime) object. The time-stamps are also converted from the local time-zone to UTC. Use data(gauges.hh.piura) to obtain an example of the format of the data-frame for which this function is designed.
#'
#' @return Object \code{gauges} A STFDF (spacetime) object containing the spatial information (i.e. gauge locations with attributes) and temporal information (i.e. precipitation time-series for each gauge location).
#'
#' @examples
#' data(gpm.hh.piura) # STFDF object with gpm imerg observations for the Piura region in Peru
#' gpm <- gpm.hh.piura
#' data(gauges.hh.piura) # STFDF object with gauge observations for the Piura region in Peru
#' gauges <- gauges.hh.piura
#' gauges <- spacetime_gauge(gauges,projection="+proj=longlat +ellps=WGS84 +datum=WGS84", timezone=-5)
#' str(gauges)
#' gauges.sp <- gauges@sp # extract spatial information of gauge locations
#' gauges.ts <- as(gauges[,,1],"xts") # extract gauges time-series

spacetime_gauge <- function(gauges,projection,timezone){

  # extraer la informacion espacial del archivo .csv y
  # convertir en un trama de datos ("data.frame")
  gauge_sp <- as.data.frame(t(gauges[1:4,2:ncol(gauges)]))
  colnames(gauge_sp) <- c("latitud","longitud","elevacion","estacion")
  gauge_sp$latitud <- as.numeric(as.character(gauge_sp$latitud))
  gauge_sp$longitud <- as.numeric(as.character(gauge_sp$longitud))
  gauge_sp$elevacion <- as.numeric(as.character(gauge_sp$elevacion))
  gauge_sp$estacion <- as.character(gauge_sp$estacion)

  ## seleccionar las coordenadas y convertir daat.frame en un objeto de mapa ("sp")
  coordinates(gauge_sp) <- c("longitud", "latitud")
  proj4string(gauge_sp) <- CRS(projection)

  ## extraer informacion temporal y convertir en un objeto de series de tiempo ("zoo" y "xts")
  gauge_z <- zoo(gauges[5:nrow(gauges),2:ncol(gauges)],
                 as.POSIXct(as.character(gauges[5:nrow(gauges),1]),tz="UTC")+(timezone*60*60)-(15*60))# adjust time zone and center time-stamp.
  gauge_df <- data.frame(as.numeric(t(gauge_z)))
  rownames(gauge_df) <- 1:nrow(gauge_df)
  colnames(gauge_df) <- "rain"

  # convert accumulation to mm/hr rates
  conversion.factor <- as.numeric(difftime(index(gauge_z)[2],index(gauge_z)[1],units="hours"))
  gauge_df <- gauge_df/conversion.factor

  # compilar los objetos sp y zoo a un objeto 3-D
  gauges <- STFDF(gauge_sp,index(gauge_z),gauge_df)

  return(gauges)

}
