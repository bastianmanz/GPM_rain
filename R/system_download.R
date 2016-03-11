#' Function to download GPM-IMERG HDF5 files from NASA http server
#'
#' @description This function downloads GPM-IMERG HDF5 files from the NASA PPS Mirador using a pre-specified list of URLs and writes these to a sub-directory of the local working directory, named ("IMERG_HDF"). Go to http://pmm.nasa.gov/data-access/downloads/gpm for more info on the products.
#'
#' @param method A character string defining which method is used for the batch download. This can be either "wget" or "curl". Users are responsible for ensuring their machine has the required functionality to run the selected method. Wget can be installed from here: http://gnuwin32.sourceforge.net/packages/wget.htm (for Windows). Ensure the executable is placed in the relevant directory.
#' @param download_list A data.frame with a single column containing character strings that define the URLs of the GPM-IMERG HDF5 files to be downloaded from the NASA Mirador http server. See data(production_url_list) for an example.
#'
#' @details This function downloads GPM-IMERG HDF5 files with global coverage as specified in the download list. Details on how to generate the download list are provided in the material accompanying this R-package. Users should check the length of the download list as the HDF5 files are about 2MB in size each. DIFFERENCE BETWEEN \code{system_download} & \code{rcurl_download}: \code{system_download} differs from \code{rcurl_download} in that it uses tools (wget, curl) external to R to download the HDF5 files from the NASA Mirador http server, whereas \code{rcurl_download} utilizes RCurl to perform the download entirely from within R. \code{rcurl_download} also obtains the data from NASA ftp servers. This requires users to register their email and supply email and password to the \code{rcurl_download} but requires no further external steps and generates its own download list within specified dates, whereas \code{system_download} requires the download list to be generated externally and loaded into R. However, no registration is required.
#'
#' @return No object is returned within R. However, the HDF5 files specified in the download list are written to a sub-directory of the local working directory, named ("IMERG_HDF").
#'
#' @examples
#' data(production_url_list) # list of GPM IMERG HDF5 files including NASA http URL
#' product <- "production"
#' method <- "wget"
#' # Not Run system_download(method,download_list=production_url_list)

system_download <- function(method="wget",download_list){

  # read file names
  filenames <- download_list
  filenames <- data.frame(filenames[1:(nrow(filenames)-1),])

  # system-based download (external to R)

  workdir <- getwd()
  dir.create(paste0(workdir,"/IMERG_HDF"))

  setwd(paste0(workdir,"/IMERG_HDF"))

  for(i in 1:nrow(filenames)){

    fi <- filenames[i,1]

    fileConn<-file("fi.txt")
    writeLines(as.character(fi), fileConn)
    close(fileConn)

    if(method=="curl"){
      system("xargs -n 1 curl -O < fi.txt") # xargs -n 1 curl -O < myfile.dat
      file.remove("fi.txt")
    } else if(method=="wget"){
      system(" wget -i fi.txt") # wget -i myfile.dat
      file.remove("fi.txt")
    }

  }

  setwd(workdir)
}


