#This version includes 0.02 subtraction from ETM+ data according to DSM's findings
## Modified on 7/05/17 by DSM

##########################################################################################
#LANDSAT PHENOLOGY ALGORITHM RUN ON EACH PIXEL WITH SUFFICIENT NO. OF OBSERVATIONS

#Phenology algorithm parameters:
#doy - Julian day of year
#yr - Year
#dYR - disturbance year 
Landsat_Phenology <- function(EVI,doy,yr,dYR){
  
  #date_txt = paste(yr,doy,sep="")
  #dup_ind = !duplicated(date_txt)
  
  #EVI = EVI[dup_ind]
  #doy = doy[dup_ind]
  #yr = yr[dup_ind]
  
  #Create "period" index for normalizing EVI time series 
  prd <- yr
  prd[yr>=1984 & yr<=1986]<-1
  prd[yr>=1987 & yr<=1989]<-2
  prd[yr>=1990 & yr<=1992]<-3
  prd[yr>=1993 & yr<=1995]<-4
  prd[yr>=1996 & yr<=1998]<-5
  prd[yr>=1999 & yr<=2001]<-6
  prd[yr>=2002 & yr<=2004]<-7
  prd[yr>=2005 & yr<=2007]<-8
  prd[yr>=2008 & yr<=2010]<-9
  prd[yr>=2011 & yr<=2014]<-10
  
  #normalized EVI threshold for long-term mean phenology date used to calculate annual phenology
  thresh <- 0.5
  
  #Initialize Phenology matrices
  phenoSPR <- matrix(0,1,31) #annual spring phenology (1984-2014)
  phenoAUT <- matrix(0,1,31) #annual autumn phenology
  ltmSPR <- matrix(0,1,1) #long-term mean spring phenology date
  ltmAUT <- matrix(0,1,1) #long-term mean autumn phenology date
  bline <- matrix(0,1,1) #long-term mean winter background EVI
  hline <- matrix(0,1,1) #long-term mean summer maximum EVI
  rsmooth <- matrix(0,1,1) #correlation between observed and smoothed EVI
  nobs <- matrix(0,1,1) #number of available EVI observations used
  
  SPR_thresh = matrix(0,1,1)
  AUT_thresh = matrix(0,1,1)
  SPRsmooth_max = matrix(0,1,1)
  SPRsmooth_min = matrix(0,1,1)
  AUTsmooth_max = matrix(0,1,1)
  AUTsmooth_min = matrix(0,1,1)
  
  info_spr <- NA
  info_aut <- NA
  
  #Process phenology dates for each pixel in EVI matrix
  
  #Apply calibration adjustment for TM vs. ETM+ (subtract 0.02 from ETM+ observations)
  #Damien's analysis
  #w <- which(sensor==7)
  #evi_row[w] <- (evi_row[w]-0.019)/1.038
  
  #Set erroneous EVI values to NA
  EVI <- as.numeric(EVI)
  EVI[EVI==-Inf | EVI==Inf]<-NA
  EVI[EVI<0 | EVI>1]<-NA
  
  #Remove NAs and, if necessary, observations after disturbance year (dYR)
  if (dYR!=0) pos <- which(is.na(EVI)==0 & yr<dYR)
  if (dYR==0) pos <- which(is.na(EVI)==0)
  EVI <- EVI[pos]
  DOY <- doy[pos]
  YR <- yr[pos]
  PRD <- prd[pos]
  nobs <- length(which(is.na(EVI)==0)) #count number of available observations
  
  if (length(EVI)>0){
    #Remove outliers using 10-day bins
    DOY_bins <- DOY
    DOY_bins <- round(DOY/10)
    b <- boxplot(EVI~DOY_bins,plot=FALSE)
    w <- which(EVI %in% b$out)
    if (length(w)>0){
      EVI <- EVI[-w]
      DOY <- DOY[-w]
      PRD <- PRD[-w]
      YR <- YR[-w]
    }
    
    
    maxEVI <- quantile(EVI,0.90)
    minEVI <- quantile(EVI,0.10)
    
    #Normalize EVI time series using 10% and 90% quantiles for EVI during three-year windows
    DF <- data.frame(PRD,EVI)
    colnames(DF) <- c('yr','evi')
    q1 <- with(DF,tapply(EVI,PRD,quantile,probs=0.10))
    q2 <- with(DF,tapply(EVI,PRD,quantile,probs=0.90))
    
    quant1 <- matrix(NA,10,1)
    quant2 <- matrix(NA,10,1)
    quant1[as.numeric(names(q1))] <- q1
    quant2[as.numeric(names(q2))] <- q2
    
    EVImax <- quant2[PRD]
    EVImin <- quant1[PRD]
    EVInorm <- (EVI-EVImin)/(EVImax-EVImin)
    
    EVInorm[EVInorm==Inf | EVInorm==-Inf]<-NA
    
    pos <- which(is.na(EVInorm)==0 & DOY<=365)
    EVInorm <- EVInorm[pos]
    DOY <- DOY[pos]
    YR <- YR[pos]
    
    if (nobs > 100){
      #Sort DOY/EVI by DOY and compute winter background EVI for improved smoothing spline fit
      x <- cbind(DOY,EVInorm) #Spring background
      x <- x[order(x[,1]),]
      start <- seq(1,x[1,1],1) 
      end <- seq(x[nrow(x),1],365,1)
      YRextra <- c(YR,rep(2012,each=length(start)),rep(2012,each=length(end)))
      DOYextra <- c(DOY,start,end)
      EVIextra <- c(EVInorm,matrix(0,length(start)+length(end),1))
      
      #Fit smoothing spline through EVI and DOY data
      fit <- smooth.spline(DOYextra,EVIextra,spar=0.55)
      EVIsmooth <- data.frame(predict(fit,x=1:365))
      rsmooth <- cor(EVIsmooth[DOY,2],EVInorm)
      
      #if (rsmooth>0.85){
      if (rsmooth>0.5){
        #Separate spline into spring and autumn segments using annual maximum      
        pkval <- which.max(EVIsmooth[,2])
        SPRsmooth <- EVIsmooth[1:pkval,2]
        AUTsmooth <- EVIsmooth[(pkval+1):365,2];
        
        #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
        #anamoly calculation
        SPR_thresh <- which.min(abs(SPRsmooth-thresh))
        SPR_halfmax <- which.min(abs(SPRsmooth-0.5))
        
        #Find anomalies inside of designated box
        SPRsmooth_max <- 1
        SPRsmooth_min <- 0
        box_max <- SPRsmooth_max-0.2*(SPRsmooth_max-SPRsmooth_min)
        box_min <- SPRsmooth_min+0.2*(SPRsmooth_max-SPRsmooth_min)
        
        #Generate a matrix with candidate spring phenology observations
        if (is.na(box_max)==0 && is.na(box_min)==0){
          info <- cbind(YR,DOY,PRD,EVInorm,matrix(0,length(DOY),1))
          info <- info[order(info[,1],info[,2]),]
          info <- rbind(matrix(0,1,5),info)
          pos <- which(is.na(info[,4])==0 & info[,4]>box_min & info[,4]<box_max &
              info[,2]<SPR_halfmax+20 & info[,2]>SPR_halfmax-20)
          if (length(pos)>4){
            k <- 1
            for (p in 1:length(pos)){
              smooth_ratio <- abs(SPRsmooth-info[pos[p],4])
              info[pos[p],5] <- which.min(smooth_ratio)
              
              if (k==1)
                info_spr <- info[pos[p],]
              else
                info_spr <- rbind(info_spr,info[pos[p],])
              
              k <- k+1
            }
            w <- which((info[pos,4]<info[(pos-1),4] & info[pos,1]==info[(pos-1),1])==1)
            if (length(w)>0) info_spr <- info_spr[-w,]
          }
        }
        
        #Compute half-maximum of spring logistic for "ruling in" image dates (points) for
        #anamoly calculation
        AUT_thresh <- which.min(abs(AUTsmooth-thresh))+pkval
        AUT_halfmax <- which.min(abs(AUTsmooth-0.5))+pkval
        
        #Find anomalies inside of designated box
        AUTsmooth_max <- 1
        AUTsmooth_min <- 0
        box_max <- AUTsmooth_max-0.3*(AUTsmooth_max-AUTsmooth_min)
        box_min <- AUTsmooth_min+0.2*(AUTsmooth_max-AUTsmooth_min)
        
        #Generate matrix with candidate autumn phenology observations
        if (is.na(box_max)==0 && is.na(box_min)==0){
          info <- cbind(YR,DOY,PRD,EVInorm,matrix(0,length(DOY),1))
          info <- info[order(info[,1],info[,2]),]
          info <- rbind(info,matrix(0,1,5))
          pos <- which(is.na(info[,4])==0 & info[,4]>box_min & info[,4]<box_max &
              info[,2]<AUT_halfmax+20 & info[,2]>AUT_halfmax-20)
          if (length(pos)>4){
            k <- 1
            for (p in 1:length(pos)){
              smooth_ratio <- abs(AUTsmooth-info[pos[p],4])
              info[pos[p],5] <- which.min(smooth_ratio)+pkval
              
              if (k==1)
                info_aut <- info[pos[p],]
              else
                info_aut <- rbind(info_aut,info[pos[p],])
              
              k <- k+1
            }
            w <- which((info[pos,4]<info[(pos+1),4] & info[pos,1]==info[(pos+1),1])==1)
            if (length(w)>0) info_aut <- info_aut[-w,]
          }
        }
        
        #Calculate interannual phenology dates by taking the distance
        #between each candidate observation and where the same magnitude 
        #of EVI occurs on the spline
        if (exists('info_spr') == 1 && exists('info_aut') == 1){
          if (length(info_spr) > 20 & length(info_aut) > 20) {
            ## dsm --- why do you use ncell here - 
            ### further why do u use ncell when you should do length(info[,1])
            ### if (ncell(info_spr) > 20 & ncell(info_aut) > 20) {
            info_spr <- cbind(info_spr,SPR_thresh+(info_spr[,2]-info_spr[,5]))
            info_aut <- cbind(info_aut,AUT_thresh+(info_aut[,2]-info_aut[,5]))
            
            for (y in 1984:2014){
              pos1 <- which(info_spr[,1] == y)
              pos2 <- which(info_aut[,1] == y)
              
              if (length(pos1) > 0 && ncol(info_spr)==6){
                w <- which.min(abs(info_spr[pos1,4]-0.5))
                phenoSPR[1,y-1983] <- ceiling(mean(info_spr[pos1[w],6]))
              }
              if (length(pos2) > 0 && ncol(info_aut)==6){
                w <- which.min(abs(info_aut[pos2,4]-0.5))
                phenoAUT[1,y-1983] <- ceiling(mean(info_aut[pos2[w],6]))
              }
            }
          }
          
          ltmSPR <- SPR_thresh
          ltmAUT <- AUT_thresh
          
          remove(info_spr,info_aut) 
        }
      }
    }
    
    pheno_matrix <- cbind(nobs,rsmooth,maxEVI,minEVI,ltmSPR,ltmAUT,phenoSPR,phenoAUT)
    
  } else {
    
    pheno_matrix <- cbind(nobs,rsmooth,NA,NA,ltmSPR,ltmAUT,phenoSPR,phenoAUT)
  }
  
#   if (is.character(pheno_matrix) == 1){
#     pheno_matrix <- matrix(NA,1,68)
#     print(paste("Pixel ID",pID,"had an error!",sep=''))
#   } 
  
  return(pheno_matrix)
}

## main function to calculate phenology metrics
run_pheno <- function(evi_dat,head_dat,dist_dat) {
  
  ## organize inputs - from header
  date_char = as.character(head_dat[,3])
  year_int = as.numeric(substring(date_char,1,4))
  date_int = as.numeric(substring(date_char,5,7))
  sens_int = head_dat[,1]
  
  all_years = seq(1984,2014)
  npix = dim(evi_dat)[[1]]
  ndates = dim(evi_dat)[[2]]
  
  ## make the lowest values/fill 0
  dist_dat[is.na(dist_dat)] = 0
  dist_dat[dist_dat<0] = 0
  
  dist_years = array(0,npix)
  dist_years[dist_dat>0] = as.numeric(substr(as.character(dist_dat[dist_dat>0]),1,4))

  pheno_out = array(NA,dim=c(npix,68))
  # pheno_out[i,] <- mapply(evi_dat, Landsat_Phenology(evi_dat[i,],date_int,year_int,dist_dat[i]) ) 
  
  ## Parallelized version
  pheno_out <- foreach(i = 1:npix, .combine = rbind) %dopar% {
    dYR <- dist_years[i] #year of first disturbance
    evi_row <- evi_dat[i,]
    pheno_matrix <- try(Landsat_Phenology(evi_row,date_int,year_int,dYR))         
  }
  
  ## Serialized version
  ##for(i in 1:100) {
     ## the inputs are the evi data, and for each pixel all the doy, year, sensor id, and disturbance value
    ##pheno_out[i,] <- try( Landsat_Phenology(evi_dat[i,],date_int,year_int,dist_dat[i]) )   
  ##}
	colnames(pheno_out) = c("nobs","rsmooth","maxEVI","minEVI","ltmSPR","ltmAUT",
				paste("spr",seq(1984,2014),sep=""),paste("aut",seq(1984,2014),sep=""))

  return(pheno_out)
} ## end run_pheno function


  
