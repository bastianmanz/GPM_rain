#' Function to perform double kriging
#'
#' @description This function performs double kriging over the satellite grid using gauge observations and climatological variograms as defined in \code{cvgm}. The process consists of two steps: (1) kriging of binary occurrence indicators and back-transformation to a binary precipitation occurrence field; and, (2) kriging of rainfall amounts using a normal-score transformation to ensure Gaussian-distributed rainfall. The binary occurrence field is used to mask the interpolated amounts field and return a combined rainfall field.
#'
#' @param gauge.ts A zoo object containing the gauge observations.
#' @param gauge.sp A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.
#' @param gpm.ts A zoo object containing the gpm-imerg observations.
#' @param gpm.sp A SpatialPointsDataFrame object containing the locations and numeric identifiers of the GPM-IMERG pixels.
#' @param ts.classification A zoo object with a single data column defining the assigned class for each time-step.
#' @param clim_vgm A list containing information on the climatological variograms as supplied by \code{cvgm}. See function \code{cvgm} for details.
#' @param rain_threshold A numeric value defining the level below which all rainfall will be ignored and measurements treated as zeroes. This level can be set to the gauge detection limit, which typically is 0.1 or 0.2 mm/hr. However, the function requires a positive value to be supplied. Hence, if only actual zero-measurements are supposed to be treated as zeroes, use a miniscule value, such as 1e-4 (default).
#' @param vgm.type A character string defining which data the variogram is to be calculated from. This can be either ordinary kriging ("OK") or kriging with external drift ("KED"; default). In the case of OK, only the gauge data is used. In the case of KED, the satellite is used as a secondary variable when defining the semi-variogram. The vgm.type should be the same as used in the fitting of the climatological semi-variogram in \code{cvgm}.
#'
#' @details In double kriging, the climatological variograms of precipitation occurrence and amounts are used to predict precipitation at pre-defined ungauged locations, which correspond to the grid nodes of the regular grid supplied by \code{gpm.sp}. In the first kriging step, binary occurrence indicators are interpolated to a continuous field ranging 0 to 1, which is then transformed into a binary field, using a threshold level of 0.5 (50% probability). In the second kriging step, amounts are interpolated and subsequently masked using the interpolated occurrence field to generate variable rainfall fields with defined spatial boundaries. For the amounts kriging, gauge observations are normal-score transformed into a Gaussian state prior to interpolation and back-transformed after interpolation and prior to masking with the binary occurrence grid. Aside from the kriging estimates, the kriging estimation variance is also calculated. For the amounts, the final kriging estimation variance is the combined (multiplied) estimation variance of the occurrence kriging and the amounts kriging. It should be noted that the back-transformation in the current version makes the kriging estimation variance not meaningful in the back-transformed state. This is to be resolved in future versions.
#'
#' @return Object \code{dk} A list containing four entries: (1) \code{dk$Zs} a time-series (zoo) object containing the estimated combined rainfall estimates (i.e. occurrence and amounts interpolation) for every grid node (as defined by \code{gpm.sp}) and time-step, (2) \code{dk$Zs_var} the combined kriging estimation variance (see details above) for each grid node and time-step, (3) \code{dk$Zof} a time-series (zoo) object containing the estimated binary rainfall occurrence field for every grid node (as defined bz \code{gpm.sp}) and time-step, (4) \code{dk$Zof_var} the binary occurrence kriging estimation variance for each grid node and timestep.
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
#' ts.classification <- classify(gauge.ts,aggregation_level="season")
#' clim.vgms.ked <- cvgm(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification,class_select=1,rain_threshold=0.001,vgm.type="KED",vgm.model=c("Exp"), class.dist=5)
#' dk.list <- DoubleKrige(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification, clim_vgm=clim.vgms.ked,rain_threshold=0.001,vgm.type="KED")
# names(dk.list)

DoubleKrige <- function(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification,clim_vgm,rain_threshold=0.001,vgm.type="KED"){

  grd <- gpm.sp
  grd_ind <- gpm.sp

  # The nscore() and backtr() functions contained in this function are the work of Ashton Shortridge and are obtained from https://www.msu.edu/~ashton/research/code/nscore.R.
  nscore <- function(x) {
    nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score
    trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
    return (list(nscore=nscore, trn.table=trn.table))
  }

  backtr <- function(scores, nscore, tails='none', draw=FALSE) {
    if(tails=='none') {   # No extrapolation
      min.x <- min(nscore$trn.table$x)
      max.x <- max(nscore$trn.table$x)
    }
    min.sc <- min(scores)
    max.sc <- max(scores)
    x <- c(min.x, nscore$trn.table$x, max.x)
    nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)
    back.xf <- approxfun(nsc,x) # Develop the back transform function
    val <- back.xf(scores)
    return(val)
  }

  Zg   <- gauge.ts
  Gdata   <- gauge.sp
  Zs   <- gpm.ts
  Tdata   <- gpm.sp

  loc <- gauge.sp@data$GPM_PIXEL

  Zg         <- data.frame(t(Zg))
  colnames(Zg) <- paste("G",1:ncol(Zg),sep="")
  gaugename  <- colnames(Zg)
  Zs_field <- as.data.frame(t(Zs))
  colnames(Zs_field) <- paste("S",1:ncol(Zs_field),sep="")
  trendname<- colnames(Zs_field)

  # Separate into Occurrence indicator and Amounts
  Zg_ind <- apply(Zg,2,function(x) replace(x,which(x>=rain_threshold),1))
  Zg_ind <- apply(Zg_ind,2,function(x) replace(x,which(x<rain_threshold),0))
  Zg <- apply(Zg,2,function(x) replace(x,which(x<rain_threshold),NA))

  Zs_field_ind <- apply(Zs_field,2,function(x) replace(x,which(x>=rain_threshold),1))
  Zs_field_ind <- apply(Zs_field_ind,2,function(x) replace(x,which(x<rain_threshold),0))
  Zs_field <- apply(Zs_field,2,function(x) replace(x,which(x<rain_threshold),NA))

  # normal-score transformation
  zg_nscore_list <- list()
  zs_nscore_list <- list()

  for(i in 1:length(gaugename)){
    if(all(is.na(Zg[,i]))==FALSE){
      test <- nscore(Zg[,i])
      Zg[,i]<- test$nscore
      zg_nscore_list[[i]]<- test
    }
    if(all(is.na(Zs_field[,i]))==FALSE){
      testzs <- nscore(Zs_field[,i])
      Zs_field[,i]<- testzs$nscore
      zs_nscore_list[[i]]<- testzs
    }
  }

  Zs_trend <- Zs_field[loc,]
  Zs_trend_ind <- Zs_field_ind[loc,]

  # amounts: event normalization
  Zg_var <- apply(Zg,2,function(x) var(x,na.rm=T))
  Zg_ind_var <- apply(Zg_ind,2,function(x) var(x,na.rm=T))
  Zs_field_var <- apply(Zs_field,2,function(x) var(x,na.rm=T))
  Zs_field_ind_var <- apply(Zs_field_ind,2,function(x) var(x,na.rm=T))
  Zs_trend_var <- apply(Zs_field,2,function(x) var(x,na.rm=T))
  Zs_trend_ind_var <- apply(Zs_field_ind,2,function(x) var(x,na.rm=T))

  if(vgm.type=="OK"){
    ind_data <- as.data.frame(Zg_ind)
    data <- as.data.frame(Zg)
    coordinates(ind_data) <- coordinates(Gdata)
    coordinates(data) <- coordinates(Gdata)

  } else if(vgm.type=="KED"){
    ind_data <- as.data.frame(cbind(Zg_ind,Zs_trend_ind))
    data <- as.data.frame(cbind(Zg,Zs_trend))
    coordinates(ind_data) <- coordinates(Gdata)
    coordinates(data) <- coordinates(Gdata)

    Zs_field <- as.data.frame(Zs_field)
    coordinates(Zs_field) <- coordinates(Tdata)
    proj4string(Zs_field) <- CRS(proj4string(Tdata))
    grd <- Zs_field
    colnames(grd@data) <- trendname

    Zs_field_ind <- as.data.frame(Zs_field_ind)
    coordinates(Zs_field_ind) <- coordinates(Tdata)
    proj4string(Zs_field_ind) <- CRS(proj4string(Tdata))
    grd_ind <- Zs_field_ind
    colnames(grd_ind@data) <- trendname

  } else if(vgm.type=="SAT"){
    ind_data <- as.data.frame(Zs_field_ind)
    data <- as.data.frame(Zs_field)
    gaugename <- trendname
    coordinates(ind_data) <- coordinates(Tdata)
    coordinates(data) <- coordinates(Tdata)
  }

  proj4string(ind_data) <- CRS(proj4string(Gdata))
  proj4string(data) <- CRS(proj4string(Gdata))

  dk <- list()
  dk$Zs          <- matrix(data=0,ncol=nrow(grd),nrow=length(gaugename))
  dk$Zs_var      <- matrix(data=0,ncol=nrow(grd),nrow=length(gaugename))
  dk$Zof         <- matrix(data=0,ncol=nrow(grd),nrow=length(gaugename))
  dk$Zof_var     <- matrix(data=0,ncol=nrow(grd),nrow=length(gaugename))

  # create vector when individual events occur
  rain_event <- apply(ind_data@data[,1:ncol(Zg)],2,function(x) sum(x==1,na.rm=T)>=3)
  full.occ.event <- apply(ind_data@data[,1:ncol(Zg)],2,function(x) sum(x==1,na.rm=T)==length(x))

  # set progress bar
  pb <- txtProgressBar(style=3)
  print("Performing double kriging interpolation over full grid")

  for(i in 1:length(rain_event)){

    if(rain_event[i]==TRUE){

      setTxtProgressBar(pb, i/max(which(rain_event==T)))

      # Get data for time step and exclude gauges with missing data
      if(vgm.type=="KED"){
        ind_data_sub    <- ind_data[,which(colnames(ind_data@data)==gaugename[i] | colnames(ind_data@data)==trendname[i])]
        data_sub    <- data[,which(colnames(data@data)==gaugename[i] | colnames(data@data)==trendname[i])]
      } else{
        ind_data_sub    <- ind_data[,which(colnames(ind_data@data)==gaugename[i])]
        data_sub    <- data[,which(colnames(data@data)==gaugename[i])]
      }

      ind_data_sub    <- ind_data_sub[which(!is.na(ind_data_sub@data[,1])),]
      data_sub    <- data_sub[which(!is.na(data_sub@data[,1])),]

      if(full.occ.event[i]==FALSE){ # i.e. do not cvgm/krige, if all gauges ==1, as computations fail
        if(vgm.type=="KED"){
          if(all(ind_data_sub@data[,2]==0) | all(ind_data_sub@data[,2]==1)){ # use OK-vgm, if all sat==0 or ==1
            ind_formula <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
          } else{ # use KED-vgm
            ind_formula <- as.formula(paste(as.character(gaugename[i])," ~ ", as.character(trendname[i]) ,sep=""))
          }
        } else if(vgm.type=="OK" | vgm.type=="SAT"){
          ind_formula <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
        }

        # define occ. event_vgms from occ. clim_vgms and scale by event var
        event_occ_vgm <- clim_vgm$occ_model[[clim_vgm$occ_event_id[i]]]
        event_occ_vgm$psill = event_occ_vgm$psill * Zg_ind_var[i]

        # Perform Indicator Kriging
        occ_fields <- krige(ind_formula,ind_data_sub,grd_ind,model=event_occ_vgm,debug.level=0)
        occ_fields@data[which(occ_fields@data[,1]<0.5),1]=0
        occ_fields@data[which(occ_fields@data[,1]>=0.5),1]=1

      } else{ # if all gauges==1, make entire occ grd ==1/ var also ==1
        occ_fields <- gpm.sp
        occ_fields@data[,1] <- 1
        occ_fields@data <- data.frame(occ_fields@data[,1],occ_fields@data[,1])
        colnames(occ_fields@data)<- c("var1.pred","var1.var")
      }

      ## AMOUNTS ##
      # define amts. svgm model
      if(vgm.type=="KED"){
        if(all(data_sub@data[,2]==0) | any(is.na(data_sub@data)) | length(unique(data_sub@data))==1){ # use OK-vgm, as all sat==0
          formula <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
        } else{ # use KED-vgm
          formula <- as.formula(paste(as.character(gaugename[i])," ~ ", as.character(trendname[i]) ,sep=""))
        }
      } else if(vgm.type=="OK" | vgm.type=="SAT"){
        formula <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
      }

      # define amts event_vgms from amts clim_vgms and scale by event var
      event_amts_vgm <- clim_vgm$amts_model[[clim_vgm$amts_event_id[i]]]
      event_amts_vgm$psill = event_amts_vgm$psill * Zg_var[i]

      # Perform Amounts Kriging
      rain_fields <- krige(formula,data_sub,grd,model=event_amts_vgm,debug.level=0)
      # back-transform
      rain_fields@data <- as.data.frame(apply(rain_fields@data,2,function(x) backtr(x,zg_nscore_list[[i]])))

      # combine occ field and amts field
      rain_fields@data[,1] <- rain_fields@data[,1]/occ_fields@data[,1]
      rain_fields@data[which(!is.finite(rain_fields@data[,1])),1]=0
      rain_fields@data[,2] <- rain_fields@data[,2]*occ_fields@data[,2]

      # store merged fields in list
      dk$Zs[i,] <- as.numeric(rain_fields@data[,1])
      dk$Zs_var[i,] <- as.numeric(rain_fields@data[,2])
      dk$Zof[i,] <- as.numeric(occ_fields@data[,1])
      dk$Zof_var[i,] <- as.numeric(occ_fields@data[,2])

    }
  }

  dk$Zs <- as.xts(zoo(dk$Zs,index(gpm.ts)))
  dk$Zs_var <- as.xts(zoo(dk$Zs_var,index(gpm.ts)))
  dk$Zof <- as.xts(zoo(dk$Zof,index(gpm.ts)))
  dk$Zof_var <- as.xts(zoo(dk$Zof_var,index(gpm.ts)))

  return(dk)

}
