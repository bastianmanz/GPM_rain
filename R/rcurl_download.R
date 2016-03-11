#' Function to download GPM-IMERG HDF5 files from NASA ftp server using RCurl
#'
#' @description This function downloads GPM-IMERG HDF5 files from the NASA PPS ftp servers and writes these to a sub-directory of the local working directory, named ("IMERG_HDF"). Unlike \code{system_download} it does nmot require a dpwnload_list to be provided externally, but creates one itself. Go to http://pmm.nasa.gov/data-access/downloads/gpm for more info on the products.
#'
#' @param product A character string defining the GPM-IMERG product to be downloaded. This can be either "production" (gauged-corrected, 4-month latency) or "nrt" (Near-real time, without gauge correction, latency < 1 day).
#' @param nrt_type A character string defining the NRT product type to be downloaded. This can be either "early" (4-hr latency, raw satellite - no correction) or "late" (18-hr latency, satellite with climatological correction). The \code{nrt_type} is only relevant if \code{product} is set to "nrt" and is ignored otherwise.
#' @param start A numeric value defining the starting date (full day) for the download in the format: YYYYMMDD (e.g. 20150401).
#' @param end A numeric value defining the ending date (full day) for the download in the format: YYYYMMDD (e.g. 20150402).
#' @param userpwd A character string defining your username and password to access the NASA ftp server. Downloading IMERG HDF5 files from the NASA ftp servers requires you to register your email here http://registration.pps.eosdis.nasa.gov/registration/newContact.html. Ensure to check the box that you are also interested in NRT products, if this is the case. Once you have done this and received the confirmation email enter your username and password in this format ("username:password"): "you(at)email:you(at)email". Your email == your username == your password.
#'
#' @details This function downloads GPM-IMERG HDF5 files from the start to the end dates from the NASA ftp servers using RCurl. DIFFERENCE BETWEEN \code{system_download} & \code{rcurl_download}: \code{system_download} differs from \code{rcurl_download} in that it uses tools (wget, curl) external to R to download the HDF5 files from the NASA Mirador http server, whereas \code{rcurl_download} utilizes RCurl to perform the download entirely from within R. \code{rcurl_download} also obtains the data from NASA ftp servers. This requires users to register their email and supply email and password to the \code{rcurl_download} but requires no further external steps and generates its own download list within specified dates, whereas \code{system_download} requires the download list to be generated externally and loaded into R. However, no registration is required.
#'
#' @return download_list A data.frame obect with a single column containing the character strings that define the file names the downloaded HDF 5 files. Additionally, the HDF5 files between the specified start and end dates are written to a sub-directory of the local working directory, named ("IMERG_HDF").
#'
#' @examples
#' data(production_url_list) # list of GPM IMERG HDF5 files including NASA http URL
#' product <- "production"
#' nrt_type <- "late"
#' start <- 20150401
#' end <- 20150402
#' userpwd <- "you(at)email:you(at)email"
#' # Not Run rcurl_download(product,nrt_type,start,end,userpwd)

rcurl_download <- function(product,nrt_type,start,end,userpwd){

  ### Obtain GPM-IMERG-NRT file names for ftp download ###

  if(product == "nrt"){

    # GENERATE DATES

    st_yr <- paste(strsplit(as.character(start),split="")[[1]][1:4],collapse="")
    st_mon <- paste(strsplit(as.character(start),split="")[[1]][5:6],collapse="")
    st_yr_mon <- as.yearmon(paste(st_yr,"-",st_mon,sep=""))
    start_day <- as.numeric(paste(strsplit(as.character(start),split="")[[1]][7:8],collapse=""))

    end_yr <- paste(strsplit(as.character(end),split="")[[1]][1:4],collapse="")
    end_mon <- paste(strsplit(as.character(end),split="")[[1]][5:6],collapse="")
    end_yr_mon <- as.yearmon(paste(end_yr,"-",end_mon,sep=""))
    end_day <- as.numeric(paste(strsplit(as.character(end),split="")[[1]][7:8],collapse=""))

    mons <- as.character(seq(from=as.Date(st_yr_mon),to=as.Date(end_yr_mon),by="month"))
    mons <- as.character(lapply(mons,function(x) paste(strsplit(x,split="")[[1]][c(1:4,6:7)],collapse="")))

    url <- paste("ftp://jsimpson.pps.eosdis.nasa.gov/data/imerg/",nrt_type,"/",mons,"/",sep="")

    ch <- getCurlHandle()

    # set progressbar
    pb <- txtProgressBar()
    print("GPM IMERG Production: Identify directories and files from ftp server...")

    for(i in 1:length(url)){
      setTxtProgressBar(pb, i/length(url))
      filenames = getURL(url[i], userpwd=userpwd,ftp.use.epsv = FALSE, dirlistonly = TRUE, curl=ch)
      dirfilenames = paste(url[i], strsplit(filenames, "\r*\n")[[1]], sep = "")
      filenames = strsplit(filenames, "\r*\n")[[1]]

      if(i==1){
        start_file <- min(which(as.numeric(lapply(filenames, function(x) paste(strsplit(strsplit(x,split="[.]")[[1]][5],split="")[[1]][7:8],collapse="")))>= start_day))
        dirfilenames <- dirfilenames[start_file:length(filenames)]
        filenames <- filenames[start_file:length(filenames)]
      }
      if(i == length(url)){
        end_file <- max(which(as.numeric(lapply(filenames, function(x) paste(strsplit(strsplit(x,split="[.]")[[1]][5],split="")[[1]][7:8],collapse="")))<= end_day))
        dirfilenames <- dirfilenames[1:end_file]
        filenames <- filenames[1:end_file]
      }
      if(exists("dirfile_list")){dirfile_list <- c(dirfile_list,dirfilenames)} else{dirfile_list <- dirfilenames}
      if(exists("file_list")){file_list <- c(file_list,filenames)} else{file_list <- filenames}
    }

    ### Obtain GPM-IMERG-PRODUCTION file names for ftp download ###

  } else if(product=="production"){

    # GENERATE DATES

    prod_dates <- seq(from=as.Date(as.character(start),format="%Y%m%d"),to=as.Date(as.character(end),format="%Y%m%d"),by="day")
    prod_dates <- as.character(format(prod_dates,"%Y/%m/%d"))

    url <- paste("ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/",prod_dates,"/imerg/",sep="")

    pb <- txtProgressBar()
    print("GPM IMERG Production: Identify directories and files from ftp server...")

    ch <- getCurlHandle()

    for(i in 1:length(url)){
      setTxtProgressBar(pb, i/length(url))
      filenames = getURL(url[i], userpwd=userpwd,ftp.use.epsv = FALSE, dirlistonly = TRUE, curl=ch)
      dirfilenames = paste(url[i], strsplit(filenames, "\r*\n")[[1]], sep = "")
      filenames = strsplit(filenames, "\r*\n")[[1]]

      if(exists("dirfile_list")){dirfile_list <- c(dirfile_list,dirfilenames)} else{dirfile_list <- dirfilenames}
      if(exists("file_list")){file_list <- c(file_list,filenames)} else{file_list <- filenames}
      rm(filenames,dirfilenames)
    }
  }

  ### BINARY DOWNLOAD AND HDF5 EXTRACTION ###

  # create sub-directory to store files in
  workdir <- getwd()
  dir.create(paste0(workdir,"/IMERG_HDF"))
  setwd(paste0(workdir,"/IMERG_HDF"))

  pb <- txtProgressBar()
  print("GPM IMERG: Binary Download from ftp server")

  ch <- getCurlHandle()

  for(i in 1:length(dirfile_list)){
    setTxtProgressBar(pb, i/length(dirfile_list))
    bin <- getBinaryURL(dirfile_list[i],userpwd=userpwd, curl=ch)
    writeBin(bin,file_list[i])
  }

  setwd(workdir) # return to org working directory

  # write download list to directory
  write.table(as.data.frame(dirfile_list),file=paste("file_list_",product,".txt",sep=""))

  # return file_list (w/o URLs)
  download_list <- as.data.frame(file_list)
  return(download_list)

}
