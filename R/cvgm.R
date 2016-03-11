#' Function to fit climatological semi-variograms
#'
#' @description This function undergoes a number of steps to calculate experimental (empirical) semi-variograms of gauge-based time-series, aggregates these based on the previously defined climatological classification approach, and, ultimately, following a normalization step, fits theoretical semi-variogram models for each class. The motivation for claculating climatological variograms lies in the fact that, depending on data availability, fitting theoretical variogram models at each individual time-step may not be successful, making climatological variograms a more robust approach. This step is essential to for Double Kriging.
#'
#' @param gauge.ts A zoo object containing the gauge observations.
#' @param gauge.sp A SpatialPointsDataFrame object containing the gauge locations and additional gauge attributes.
#' @param gpm.ts A zoo object containing the gpm-imerg observations.
#' @param gpm.sp A SpatialPointsDataFrame object containing the locations and numeric identifiers of the GPM-IMERG pixels.
#' @param ts.classification A zoo object with a single data column defining the assigned class for each time-step.
#' @param class_select A positive integer defining the column from the object \code{ts.classification} from which the climatological classification is to be obtained.
#' @param rain_threshold A numeric value defining the level below which all rainfall will be ignored and measurements treated as zeroes. This level can be set to the gauge detection limit, which typically is 0.1 or 0.2 mm/hr. However, the function requires a positive value to be supplied. Hence, if only actual zero-measurements are supposed to be treated as zeroes, use a miniscule value, such as 1e-4 (default).
#' @param vgm.type A character string defining which data the variogram is to be calculated from. This can be either ordinary kriging (OK), kriging with external drift (KED; default) or the satellite (SAT). In the case of OK, only the gauge data is used. In the case of KED, the satellite is used as a secondary variable when defining the semi-variogram. In the case of SAT the semi-variogram is calculated from the full grid of satellite observations. If subsequently using the semi-variograms for kriging, use the same vgm.type in cvgm and in DK.
#' @param vgm.model A character string defining the theoretical semi-variogram model that is to be fitted to the averaged and normalized experimental semi-variogram. These models are derived from the gstat package. Type vgm() to obtain a full list. Multiple models can be provided and the one with the lowest sum of squared errors is selected and returned.
#' @param class.dist A numeric value defining the width (horizontal distance) in kilometres over which the experimental variograms within a climatological class are binned. Increasing this value will smooth out the variability and improve the fitting of the theoretical variogram, but may lead to numerical artefacts. A low value typically results in a lot of noise. The value should not be larger than the grid cell size fo the satellite product/ target grid (0.1 degree in the case of GPM-IMERG).
#'
#' @details The climatological variogram function is a required pre-requisite for the Double Kriging (DK) step as it defines the spatial variance structure (semi-variogram) that is required for prediction in DK. The input time-series data is converted into an occurrence time-series (binary indicator) or a normal-score transformed Gaussian rainfall intensity time-series. For both of these time-series an experimental semi-variogram is determined at each time-step. After this is done for the full time-series, the semi-variograms are classified depending on the pre-defined climatological classification and normalized and averaged for each class.
#'
#' @return Object \code{clim_vgm} A list containing four entries for each, occurrence and amounts. An event identifier is provided that defines the climatological class for event time-step/ event. The second entry contains the fitted models required for the kriging. The third entry is a data frame summarizing the climatological semi-variogram and number of events for each climatological class. The last entry contains the semi-variogram plots. Type names(clim_vgm) to obtain the names of all entries.
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
#' clim.vgms.ked <- cvgm(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification,class_select=1,rain_threshold=0.001,vgm.type="KED",vgm.model=c("Exp","Gau","Sph","Mat"), class.dist=5)
#' class(clim.vgms.ked)
#' names(clim.vgms.ked)
#' clim.vgms.ked$occ_vgm_results
#' clim.vgms.ked$amts_vgm_results
#' print(clim.vgms.ked$occ_plot[[1]], split=c(1,1,4,2), more=TRUE)
#' print(clim.vgms.ked$occ_plot[[2]], split=c(2,1,4,2), more=TRUE)
#' print(clim.vgms.ked$occ_plot[[3]], split=c(3,1,4,2), more=TRUE)
#' print(clim.vgms.ked$occ_plot[[4]], split=c(4,1,4,2), more=TRUE)
#' print(clim.vgms.ked$amts_plot[[1]], split=c(1,2,4,2), more=TRUE)
#' print(clim.vgms.ked$amts_plot[[2]], split=c(2,2,4,2), more=TRUE)
#' print(clim.vgms.ked$amts_plot[[3]], split=c(3,2,4,2), more=TRUE)
#' print(clim.vgms.ked$amts_plot[[4]], split=c(4,2,4,2), more=TRUE)

cvgm <- function(gauge.ts,gauge.sp,gpm.ts,gpm.sp,ts.classification,class_select=1,rain_threshold=0.0001,vgm.type="KED",vgm.model="Exp", class.dist=5){

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
  Tdata 	<- gpm.sp

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

  } else if(vgm.type=="SAT"){
    ind_data <- as.data.frame(Zs_field_ind)
    data <- as.data.frame(Zs_field)
    gaugename <- trendname
    coordinates(ind_data) <- coordinates(Tdata)
    coordinates(data) <- coordinates(Tdata)
  }

  proj4string(ind_data) <- CRS(proj4string(Gdata))
  proj4string(data) <- CRS(proj4string(Gdata))

  emp_vgm_occ <- list()
  emp_vgm_amts <- list()

  # create vector when individual events occur
  rain_event <- apply(ind_data@data[,1:ncol(Zg)],2,function(x) sum(x==1,na.rm=T)>=3)
  full.occ.event <- apply(ind_data@data[,1:ncol(Zg)],2,function(x) sum(x==1,na.rm=T)==length(x))

  # set progress bar
  pb <- txtProgressBar(style=3)
  print("Calculating experimental semi-variograms for each event")

  for (i in as.numeric(which(rain_event==T))){
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

    ## OCCURRENCE INDICATOR ##
    # define occ.ind. emp. svgm model eqn
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

      # Model occ. ind. semivariogram
      if(exists("ind_formula")){ # i.e. if vgm model has been defined
        emp_vgm_occ[[i]] <- variogram(ind_formula, data=ind_data_sub)
        if(length(variogram(ind_formula, data=ind_data_sub))>0){ # i.e. if emp vgm. has been determined
          emp_vgm_occ[[i]]$gamma <- emp_vgm_occ[[i]]$gamma/Zg_ind_var[i]# clim_vgm = vgm/var
          emp_vgm_occ[[i]]$dist <- round(emp_vgm_occ[[i]]$dist/class.dist,digits=0)*class.dist # standardise individual emp vgms distances
        }
      }
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

    # Model amts. semivariogram
    if(exists("formula")){ # i.e. if vgm model has been defined
      emp_vgm_amts[[i]] <- variogram(formula, data=data_sub)
      if(length(variogram(formula, data=data_sub))>0){ # i.e. if emp vgm. has been determined
        emp_vgm_amts[[i]]$gamma <- emp_vgm_amts[[i]]$gamma/Zg_var[i]
        emp_vgm_amts[[i]]$dist <- round(emp_vgm_amts[[i]]$dist/class.dist,digits=0)*class.dist
      }
    }

    rm(ind_formula,formula,ind_data_sub,data_sub)
  }

  ######################################################
  # define classes for clim vgms of occurrence indicator

  classification <- ts.classification
  classes <- data.frame(sort(unique(classification)))
  colnames(classes) <- colnames(classification)

  event_classes <- list()
  class_emp_vgm <- list()

  vgm_models <- vgm.model
  fitted_vgms <- list()
  fvgm_sserr <- rep(NA,2*length(vgm_models))
  vm.fit_occ    <- list()

  # set progress bar
  pb <- txtProgressBar()
  print("Fitting occurence semi-variogram model for each climatological class")

  # compile and weighted avg empirical vgms of occ ind for each class
  for(i in 1:nrow(classes)){
    setTxtProgressBar(pb, i/nrow(classes))
    event_classes[[i]] <- as.numeric(which(rain_event==T))[as.numeric(which(apply(as.data.frame(classification[as.numeric(which(rain_event==T)),class_select]),1,function(x) all(as.numeric(x)==as.numeric(classes[i,])))==TRUE))]
    if(length(event_classes[[i]])>1){
      emp_vgms_compiled <- do.call(rbind,emp_vgm_occ[event_classes[[i]]])
      emp_vgms_compiled <- emp_vgms_compiled[with(emp_vgms_compiled,order(dist)),]
      emp_vgms_compiled <- emp_vgms_compiled[which(emp_vgms_compiled$gamma<2),]
      for(j in 1:length(unique(emp_vgms_compiled$dist))){
        if(j==1){
          class_emp_vgm[[i]] <- emp_vgms_compiled[1,]
          id_j <- which(emp_vgms_compiled$dist==unique(emp_vgms_compiled$dist)[j])
          class_emp_vgm[[i]]$np[j] <- round(mean(emp_vgms_compiled$np[id_j]))
          class_emp_vgm[[i]]$dist[j] <- unique(emp_vgms_compiled$dist)[j]
          class_emp_vgm[[i]]$gamma[j] <- sum(emp_vgms_compiled$gamma[id_j]*(emp_vgms_compiled$np[id_j]/sum(emp_vgms_compiled$np[id_j])))
        } else{
          class_emp_vgm[[i]] <- rbind(class_emp_vgm[[i]],emp_vgms_compiled[1,])
          id_j <- which(emp_vgms_compiled$dist==unique(emp_vgms_compiled$dist)[j])
          class_emp_vgm[[i]]$np[j] <- round(mean(emp_vgms_compiled$np[id_j]))
          class_emp_vgm[[i]]$dist[j] <- unique(emp_vgms_compiled$dist)[j]
          class_emp_vgm[[i]]$gamma[j] <- sum(emp_vgms_compiled$gamma[id_j]*(emp_vgms_compiled$np[id_j]/sum(emp_vgms_compiled$np[id_j])))
        }
      }
      # find opt fit for each class emp vgm
      for(k in 1:length(vgm_models)){
        fitted_vgms[[k]] <- fit.variogram(class_emp_vgm[[i]],model=vgm(psill=mean(c(max(class_emp_vgm[[i]]$gamma),median(class_emp_vgm[[i]]$gamma))),model=vgm_models[k],range=30,nugget=min(class_emp_vgm[[i]]$gamma)))
        fvgm_sserr[k] <- attr(fitted_vgms[[k]],"SSErr")
        fitted_vgms[[k+length(vgm_models)]] <- fit.variogram(class_emp_vgm[[i]],model=vgm(psill=mean(c(max(class_emp_vgm[[i]]$gamma),median(class_emp_vgm[[i]]$gamma))),model=vgm_models[k],range=30))
        fvgm_sserr[k+length(vgm_models)] <- attr(fitted_vgms[[k+length(vgm_models)]],"SSErr")
      }
      vm.fit_occ[[i]] <- fitted_vgms[[which.min(fvgm_sserr)]]
    }
  }

  num_events <- unlist(lapply(event_classes,function(x) length(x)))
  num_events <- as.data.frame(cbind(classes,num_events))

  vgm_results <- as.data.frame(matrix(data=NA,nrow=nrow(classes),ncol=4))
  colnames(vgm_results) <- c("nugget","sill","range","sserr")
  for(i in 1:length(vm.fit_occ)){
    vgm_results[i,] <- c(vm.fit_occ[[i]]$psill[1],sum(vm.fit_occ[[i]]$psill[1:2]),vm.fit_occ[[i]]$range[2],attr(vm.fit_occ[[i]],"SSErr"))
  }
  vgm_results <- cbind(num_events,vgm_results)

  event_id <- rep(NA,length(gaugename))
  for(i in 1:length(event_classes)){
    event_id[event_classes[[i]]]=i
  }

  clim_vgm <- list()
  clim_vgm$occ_event_id <- event_id
  clim_vgm$occ_model <- vm.fit_occ
  clim_vgm$occ_vgm_results <- vgm_results
  clim_vgm$occ_plot <- list()

  sns <- c("DJF","MAM","JJA","SON")

  for(i in 1:nrow(classes)){
    clim_vgm$occ_plot[[i]] <- plot(class_emp_vgm[[i]],
                                   vm.fit_occ[[i]],
                                   ylim=c(0,2),
                                   main=paste("Occ ",
                                              vgm.type," ",
                                              sns[num_events$season[i]]," ",
                                              as.character(clim_vgm$occ_model[[1]][2,1]),
                                              ", n=",num_events$num_events[i],
                                              ", SSE=",round(clim_vgm$occ_vgm_results$sserr[i],digits=3),sep=""))
  }

  ##############################################
  # define classes for clim vgms of rainfall amounts

  event_classes <- list()
  class_emp_vgm <- list()

  vgm_models <- vgm.model
  fitted_vgms <- list()
  fvgm_sserr <- rep(NA,2*length(vgm_models))
  vm.fit_amts    <- list()

  # set progress bar
  pb <- txtProgressBar()
  print("Fitting rainfall amounts semi-variogram model for each climatological class")

  # compile and weighted avg empirical vgms of rain amts for each class
  for(i in 1:nrow(classes)){
    setTxtProgressBar(pb, i/nrow(classes))
    event_classes[[i]] <- as.numeric(which(rain_event==T))[as.numeric(which(apply(as.data.frame(classification[as.numeric(which(rain_event==T)),class_select]),1,function(x) all(as.numeric(x)==as.numeric(classes[i,])))==TRUE))]
    if(length(event_classes[[i]])>1){
      emp_vgms_compiled <- do.call(rbind,emp_vgm_amts[event_classes[[i]]])
      emp_vgms_compiled <- emp_vgms_compiled[with(emp_vgms_compiled,order(dist)),]
      emp_vgms_compiled <- emp_vgms_compiled[which(emp_vgms_compiled$gamma<2),]
      for(j in 1:length(unique(emp_vgms_compiled$dist))){
        if(j==1){
          class_emp_vgm[[i]] <- emp_vgms_compiled[1,]
          id_j <- which(emp_vgms_compiled$dist==unique(emp_vgms_compiled$dist)[j])
          class_emp_vgm[[i]]$np[j] <- round(mean(emp_vgms_compiled$np[id_j]))
          class_emp_vgm[[i]]$dist[j] <- unique(emp_vgms_compiled$dist)[j]
          class_emp_vgm[[i]]$gamma[j] <- sum(emp_vgms_compiled$gamma[id_j]*(emp_vgms_compiled$np[id_j]/sum(emp_vgms_compiled$np[id_j])))
        } else{
          class_emp_vgm[[i]] <- rbind(class_emp_vgm[[i]],emp_vgms_compiled[1,])
          id_j <- which(emp_vgms_compiled$dist==unique(emp_vgms_compiled$dist)[j])
          class_emp_vgm[[i]]$np[j] <- round(mean(emp_vgms_compiled$np[id_j]))
          class_emp_vgm[[i]]$dist[j] <- unique(emp_vgms_compiled$dist)[j]
          class_emp_vgm[[i]]$gamma[j] <- sum(emp_vgms_compiled$gamma[id_j]*(emp_vgms_compiled$np[id_j]/sum(emp_vgms_compiled$np[id_j])))
        }
      }
      # find opt fit for each class emp vgm
      for(k in 1:length(vgm_models)){
        if(min(class_emp_vgm[[i]]$dist)==0){
          class_emp_vgm[[i]] <- class_emp_vgm[[i]][2:nrow(class_emp_vgm[[i]]),]
        }
        fitted_vgms[[k]] <- fit.variogram(class_emp_vgm[[i]],model=vgm(psill=mean(c(max(class_emp_vgm[[i]]$gamma),median(class_emp_vgm[[i]]$gamma))),model=vgm_models[k],range=50,nugget=min(class_emp_vgm[[i]]$gamma)))
        fvgm_sserr[k] <- attr(fitted_vgms[[k]],"SSErr")
        fitted_vgms[[k+length(vgm_models)]] <- fit.variogram(class_emp_vgm[[i]],model=vgm(psill=mean(c(max(class_emp_vgm[[i]]$gamma),median(class_emp_vgm[[i]]$gamma))),model=vgm_models[k],range=50))
        fvgm_sserr[k+length(vgm_models)] <- attr(fitted_vgms[[k+length(vgm_models)]],"SSErr")
      }
      vm.fit_amts[[i]] <- fitted_vgms[[which.min(fvgm_sserr)]]
    }
  }

  num_events <- unlist(lapply(event_classes,function(x) length(x)))
  num_events <- as.data.frame(cbind(classes,num_events))

  vgm_results <- as.data.frame(matrix(data=NA,nrow=nrow(classes),ncol=4))
  colnames(vgm_results) <- c("nugget","sill","range","sserr")
  for(i in 1:length(vm.fit_occ)){
    vgm_results[i,] <- c(vm.fit_amts[[i]]$psill[1],sum(vm.fit_amts[[i]]$psill[1:2]),vm.fit_amts[[i]]$range[2],attr(vm.fit_amts[[i]],"SSErr"))
  }
  vgm_results <- cbind(num_events,vgm_results)

  event_id <- rep(NA,length(gaugename))
  for(i in 1:length(event_classes)){
    event_id[event_classes[[i]]]=i
  }

  clim_vgm$amts_event_id <- event_id
  clim_vgm$amts_model <- vm.fit_amts
  clim_vgm$amts_vgm_results <- vgm_results
  clim_vgm$amts_plot <- list()

  for(i in 1:nrow(classes)){
    clim_vgm$amts_plot[[i]] <- plot(class_emp_vgm[[i]],
                                    vm.fit_amts[[i]],
                                    ylim=c(0,2),
                                    main=paste("Amts ",
                                               vgm.type," ",
                                               sns[num_events$season[i]]," ",
                                               as.character(clim_vgm$amts_model[[1]][2,1]),
                                               ", n=",num_events$num_events[i],
                                               ", SSE=",round(clim_vgm$amts_vgm_results$sserr[i],digits=3),sep=""))
  }

  return(clim_vgm)
}
