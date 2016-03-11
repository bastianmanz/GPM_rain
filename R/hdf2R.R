#' Function to convert GPM-IMERG files from .HDF5 to .RData
#'
#' @description This function reads GPM IMERG .HDF5 files and extracts the variables "precipitationCal", "lat" and "lon" from all the specified HDF5 files and converts these into a STFDF (spacetime) object in R, which is written to the local directory as an .RData file and also returned. See NASA GPM Technical Documentation for format and structure of HDF5 files and meaning of variables. Also, see lines 59-75 below on current structure (GPM-IMERG v3, Feb. 2016). Be aware of the size of HDF5 when specifying the spatial bbox: as GPM IMERG is half-hourly, 48 values are stored for every grid node for every day. This can results in large amounts of data being stored in memory when working at national scales.
#'
#' @param product A character string defining the GPM-IMERG product to be downloaded. This can be either "production" (gauged-corrected, 4-month latency) or "nrt" (Near-real time, without gauge correction, latency < 1 day).
#' @param nrt_type A character string defining the NRT product type to be downloaded. This can be either "early" (4-hr latency, raw satellite - no correction) or "late" (18-hr latency, satellite with climatological correction). The \code{nrt_type} is only relevant if \code{product} is set to "nrt" and is ignored otherwise.
#' @param bbox A numeric vector containing four values that define the limits of the spatial bounding box in degrees latitude/longitude. The values are to be netered in this order: minLat, maxLat, minLon, maxLon.
#' @param download_type A character string that defines which method was used to download the HDF5 files (i.e. \code{system_download} or \code{rcurl_download}). This affects the character strings in the download_list.txt file.
#'
#' @details This function reads GPM IMERG .HDF5 files and extracts the variables "precipitationCal", "lat" and "lon" from all the specified HDF5 files and converts these into a STFDF (spacetime) object in R, which is written to the local directory as an .RData file and also returned. Be aware of the size of HDF5 when specifying the spatial bbox: as GPM IMERG is half-hourly, 48 values are stored for every grid node for every day. This can results in large amounts of data being stored in memory when working at national scales.
#'
#' @return gpm A STFDF (spacetime) object containing the spatial information (grid node locations) and temporal information (observational time-series) of the extracted GPM-IMERG product. The object is also written to the local directory in an .RData file.
#'
#' @examples
#' data(production_url_list) # list of GPM IMERG HDF5 files including NASA http URL
#' product <- "production"
#' nrt_type <- "late"
#' bbox <- c(-6,-4,-80,-78)
#' download_type <-"system"
#' # Not Run gpm <- hdf2R(product, nrt_product, bbox,download_type)

hdf2R <- function(product,nrt_type,bbox,download_type){

  print(paste("Checking for file_list_",product,".txt ...",sep=""))

  # prepare .hdf5 file names for reading into R
  file_list <- read.table(paste("file_list_",product,".txt",sep=""))

  if(download_type=="system"){
    file_names <- strsplit(as.character(file_list[1,1]),split="/")[[1]][11]
    for(i in 2:nrow(file_list)){
      file_name_i <- strsplit(as.character(file_list[i,1]),split="/")[[1]][11]
      if(length(strsplit(file_name_i,split="")[[1]]) == 60){file_names <- c(file_names,file_name_i)}
    }
    file_names <- file_names[rev(1:length(file_names))]

  } else if(download_type=="rcurl"){
    file_names <- strsplit(as.character(file_list[1,1]),split="/")[[1]][9]
    for(i in 2:nrow(file_list)){
      file_name_i <- strsplit(as.character(file_list[i,1]),split="/")[[1]][9]
      if(length(strsplit(file_name_i,split="")[[1]]) == 60){file_names <- c(file_names,file_name_i)}
    }
  }

  # go to IMERG_HDF sub-directory to read files
  workdir <- getwd()
  setwd(paste0(workdir,"/IMERG_HDF"))

  south <- bbox[1]
  north <- bbox[2]
  west <- bbox[3]
  east <- bbox[4]

  imerg_cal_list <- list()

  # set progress bar
  pb <- txtProgressBar(style=3)
  print("Opening .HDF5 files and extracting spatial domain.")

  for(i in 1:length(file_names)){

    setTxtProgressBar(pb, i/length(file_names))

    # View Structure of HDF5 file #
    # h5ls(file_names[i]) # uncomment this line
    #
    # FILE STRUCTURE (DATA FIELDS)
    # group                           name       otype  dclass         dim
    # 0      /                           Grid   H5I_GROUP
    # 1  /Grid              HQobservationTime H5I_DATASET INTEGER 1800 x 3600
    # 2  /Grid                 HQprecipSource H5I_DATASET INTEGER 1800 x 3600
    # 3  /Grid                HQprecipitation H5I_DATASET   FLOAT 1800 x 3600
    # 4  /Grid           IRkalmanFilterWeight H5I_DATASET INTEGER 1800 x 3600
    # 5  /Grid                IRprecipitation H5I_DATASET   FLOAT 1800 x 3600
    # 6  /Grid                            lat H5I_DATASET   FLOAT        1800
    # 7  /Grid                            lon H5I_DATASET   FLOAT        3600
    # 8  /Grid               precipitationCal H5I_DATASET   FLOAT 1800 x 3600
    # 9  /Grid             precipitationUncal H5I_DATASET   FLOAT 1800 x 3600
    # 10 /Grid probabilityLiquidPrecipitation H5I_DATASET INTEGER 1800 x 3600
    # 11 /Grid                    randomError H5I_DATASET   FLOAT 1800 x 3600

    # read hdf5 file into R
    cal_prec <- h5read(paste(getwd(),"/",file_names[i],sep=""),name="/Grid/precipitationCal")
    lat <- h5read(paste(getwd(),"/",file_names[i],sep=""),name="/Grid/lat")
    lon <- h5read(paste(getwd(),"/",file_names[i],sep=""),name="/Grid/lon")

    # subset pixels
    pixels_lat <- which(lat>=south & lat<=north)
    pixels_lon <- which(lon>=west & lon<=east)
    imerg_cal <- cal_prec[pixels_lat,pixels_lon]
    imerg_cal <- as.numeric(t(imerg_cal))

    # time stamp
    date <- strsplit(as.character(file_names[i]),split="[.]")[[1]][5]
    date <- strsplit(date,split="-")[[1]][1]
    date <- strsplit(date,split="")[[1]]
    date <- paste(date[1],date[2],date[3],date[4],"-",date[5],date[6],"-",date[7],date[8],sep="")

    if(product=="production"){
      time <- strsplit(as.character(file_names[i]),split="-")[[1]][3]
    } else if(product=="nrt"){
      time <- strsplit(as.character(file_names[i]),split="-")[[1]][4]
    }
    time <- strsplit(time,split="")[[1]]
    time <- paste(time[2],time[3],":",time[4],time[5],":",time[6],time[7],sep="")

    date_time <- as.POSIXct(paste(date,time),tz="UTC")+(15*60) # add 15min to centre time-stamp on obs. window

    # build time-series
    if(exists("time_series")){time_series <- c(time_series,date_time)}else{time_series <- date_time}
    if(exists("imerg_cal_ts")){imerg_cal_ts <- c(imerg_cal_ts,i)}else{imerg_cal_ts <- i}
    imerg_cal_list[[i]] <- imerg_cal

    gc()
  }

  setwd(workdir) # return to org working directory

  # ensure dates and data are ordered chronologically
  names(imerg_cal_list) <- imerg_cal_ts
  df <- data.frame(time_series,imerg_cal_ts)
  df <- df[with(df, order(time_series)), ]
  imerg_cal_ts <- unlist(imerg_cal_list[df$imerg_cal_ts])
  time_series <- df$time_series

  lat <- seq(-89.95,89.95,by=0.1)
  imerg_lat <- lat[which(lat>=south & lat<=north)]
  lon <- seq(-179.95,179.95,by=0.1)
  imerg_lon <- lon[which(lon>=west & lon<=east)]
  imerg_grd <- expand.grid(x = imerg_lon, y = imerg_lat)
  colnames(imerg_grd)<-c("lon","lat")
  imerg_grd <- SpatialPoints(imerg_grd)
  imerg_grd = SpatialPointsDataFrame(imerg_grd, data.frame(ID=c(1:length(imerg_grd))))
  proj4string(imerg_grd) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

  imerg_cal_ts[imerg_cal_ts<0]=NA
  imerg_cal_ts <- data.frame(imerg_cal_ts)
  gpm <- STFDF(imerg_grd,time_series,imerg_cal_ts) # build STFDF
  print(paste("Writing file gpm_imerg_",product,".RData",sep=""))
  save(gpm,file=paste("gpm_imerg_",product,".RData",sep=""))
  rm(time_series,imerg_cal_ts)
  gc()

  return(gpm)

}
