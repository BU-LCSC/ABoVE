## main function to load the netcdf file for each row chunk
load_nc_file <- function(cur_file) {
  nc_dat = nc_open(cur_file)
  #print(nc_dat)

  head_dat = ncvar_get(nc_dat,"head")
  ts_dat = ncvar_get(nc_dat,"data")
  nc_close(nc_dat)

  nbands=dim(ts_dat)[[1]]
  ndates=dim(ts_dat)[[2]]
  ncols=dim(ts_dat)[[3]]
  nrows=dim(ts_dat)[[4]]
  npix = ncols*nrows

  all_dates = head_dat[3,]
  if(length(all_dates)!=ndates) {
    print(paste("Problem with data.
                Should be ",length(all_dates)," files but only ",ndates,".",sep=""))
  }

  ## sort the dates and save the order
  ## so that it can be used to rearrange the reflectance data
  date_order = order(all_dates)
  head_dat = t(head_dat[,date_order])
  all_dates = head_dat[,3]

  # choose the bands we want - blue, red, nir, swir1, swir2
  extract.bands <- c(1,3,4,5,6)
  num_bands = length(extract.bands)

  out_list = vector("list",num_bands)
  band_count = 0
  for(b in 1:num_bands) {
    out_list[[b]] = array(NA,dim=c(npix,ndates))
    for(y in 1:nrows) {
      start = ((y-1)*ncols) + 1
      end = y*ncols
      ## we transpose them so that pixels are rows and cols are dates
      out_list[[b]][c(start:end),] = t(ts_dat[extract.bands[b],date_order,,y]/10000)
    }
    ## add the time headers to each file
    colnames(out_list[[b]]) = all_dates
  }
  ## get fmask data too
  fmask = array(NA,dim=c(npix,ndates))
  for(y in 1:nrows) {
    start = ((y-1)*ncols) + 1
    end = y*ncols
    ## we transpose them so that pixels are rows and cols are dates
    fmask[c(start:end),] = t(ts_dat[8,date_order,,y])
  }
  rm(ts_dat)

  names(out_list) = paste("b",extract.bands,sep="")

 ## here we replace fill values (-9999 or fmask > 1) with NAs
  na_ind= fmask>1
  for(b in 1:num_bands) {
    ##
    band_bad = out_list[[b]]<0 | out_list[[b]]>1
    out_list[[b]][na_ind & band_bad]=NA
  }
## now it is ready to go - each item in the list is a band (1,2,3,4,5,7,6,fmask - or whichever you want)
## each col in the output is a single date and each row is a pixel
## the cols should correspond to the dates in head_dat - sorted
## to see the time series for a pixel - reference it as out_list[[band]][pix,]

  return(list(out_list,head_dat))
}   ## end load nc file function

## main function to composite peak summer ndvi
comp_ndvi <- function(ref_list,head_dat) {

    etm_ind = head_dat[,1]==7
    tm5_ind = head_dat[,1]==5

    full_date = head_dat[,3]
    date_char = as.character(full_date)
    year_int = as.numeric(substring(date_char,1,4))
    date_int = as.numeric(substring(date_char,5,7))

    date_ind = date_int > 179 & date_int < 241
    all_years_tm5 = seq(1984,2011)
    all_years_etm = seq(1999,2014)
    num_years_tm5 = length(all_years_tm5)
    num_years_etm = length(all_years_etm)
    npix = dim(ref_list[["red"]])[[1]]

    ## setup output list - will have one etm and one tm matrix with different dimensions
    green_out = vector("list",2)

    names(green_out) = c("etm","tm5")

    ## only interested in peak summer obs
    years_etm = year_int[date_ind & etm_ind]
    dates_etm = date_int[date_ind & etm_ind]

    years_tm5 = year_int[date_ind & tm5_ind]
    dates_tm5 = date_int[date_ind & tm5_ind]

    num_dates_tm5 = length(years_tm5)
    num_dates_etm = length(years_etm)

    band_names = c("red","nir","swir1","swir2")

    green_out[["tm5"]] <- foreach(i = 1:npix, .combine = rbind) %dopar% {
    	#print(i)

        ## make a matrix that is ndates rows by nbands cols
        cur_ref = cbind(ref_list[["red"]][i,date_ind & tm5_ind],ref_list[["nir"]][i,date_ind & tm5_ind],
                        ref_list[["swir1"]][i,date_ind & tm5_ind],ref_list[["swir2"]][i,date_ind & tm5_ind])
        colnames(cur_ref) = band_names
        cur_bool = red_filter(cur_ref[,"red"])

        cur_ref[!cur_bool,] = NA

        max_ndvi_num_yr(cur_ref,years_tm5,dates_tm5,all_years_tm5,1)
    }

    green_out[["etm"]] <- foreach(i = 1:npix, .combine = rbind) %dopar% {
      cur_ref = cbind(ref_list[["red"]][i,date_ind & etm_ind],ref_list[["nir"]][i,date_ind & etm_ind],
                      ref_list[["swir1"]][i,date_ind & etm_ind],ref_list[["swir2"]][i,date_ind & etm_ind])
      colnames(cur_ref) = band_names
      cur_bool = red_filter(cur_ref[,"red"])
      cur_ref[!cur_bool,]=NA
      ## last argument is the number of years to composite
      max_ndvi_num_yr(cur_ref,years_etm,dates_etm,all_years_etm,1)
    }

    return(green_out)

}

## composites_ndvi for each year or for a variable period of time
max_ndvi_num_yr <- function(ref_dat,years,dates,all_years,num_int=1) {

    start_year = min(all_years)
    end_year = max(all_years)
    out_years = seq(start_year,end_year,num_int)
    num_bands = dim(ref_dat)[[2]]
    num_years = length(all_years)
    ## we add two for the qa and the date
    out = array(NA,(num_years*(num_bands+2)))
    for(y in 1:num_years)
    {
        ## select only the data for that year that isnt NA
        temp_ind = years >= out_years[y] & years < out_years[y]+num_int
        na_ind = !is.na(ref_dat[,"red"])
        cur_ref = ref_dat[temp_ind & na_ind,]
        cur_dates = dates[temp_ind & na_ind]
        num_dates = length(cur_dates)
        if(num_dates > 0)
        {
          ## remove duplicated values - happens sometimes with overlapping paths
          #no_dup_ind = !duplicated(cur_dates)
          #cur_dates = cur_dates[no_dup_ind]
          out_start = ((y-1)*(num_bands+2)) + 1
          out_end = y*(num_bands+2)

          if(length(cur_dates)>1) {
            #cur_ref = cur_ref[no_dup_ind,]
            ndvi_dat <- (cur_ref[,"nir"]-cur_ref[,"red"])/(cur_ref[,"nir"]+cur_ref[,"red"])
            ndvi_dat = ndvi_dat[(ndvi_dat>0 & ndvi_dat<1)]
            ## get max of what is left
            cur_ind = get_which_max(ndvi_dat)
            out[out_start:out_end] = c(cur_ref[cur_ind,],num_dates,cur_dates[cur_ind])
          }
        } else {
          if(num_dates == 1) {
            out[out_start:out_end] = c(cur_ref,num_dates,cur_dates)
          }
        }
    }

    return(out)
}

# Hampel median filter with a custom min/max threshold
## for filtering outliers - need library(pracma)
hampel_cust <- function (x, k, t0 = 3) {
  n <- length(x)
  y <- x
  ind <- c()
  L <- 1.4826
  for (i in (k + 1):(n - k)) {
    x0 <- median(x[(i - k):(i + k)])
    S0 <- L * median(abs(x[(i - k):(i + k)] - x0))
    if ((x[i] - x0) > (t0 * S0)) {
      y[i] <- x0
      ind <- c(ind, i)
      ## added 6/29/16 bigger threshold on low outliers
    } else if((x[i] - x0) < (-1*(t0+0.5) * S0)) {
      y[i] <- x0
      ind <- c(ind, i)
    }
  }
  list(y = y, ind = ind)
}

## now we have all the red values and we produce a boolean array with the good values at each pixel
red_filter <- function(red_vals) {

  ## make indices to remove the duplicated values and the values with NAs
  #dup_ind = !duplicated(names(red_vals))
  na_ind = !is.na(red_vals)

  good_ind = red_vals>0 & red_vals<1

  ## this will be the output array
  red_bool = array(F,length(red_vals))

  ## need to remove NA data for filter to work
  #all_good = na_ind & dup_ind & good_ind
  all_good = na_ind & good_ind
  temp_ts = red_vals[all_good]

  act_len = length(temp_ts)

  orig_ts = temp_ts

  ## buffer will be the size of the buffer on either side as well as the size of the window
  ## make longer buffer - 7/3/17
  buffer = 16
  ## only do the filter if we have at least 3 good values
  if(act_len > 2) {
    #temp_med = median(temp_ts)
    ## actually use the 0.33 quantile
    temp_med = quantile(temp_ts,0.33)
    ## we put fill values on the side so the median filter will work properly
    temp_ts = c(rep(temp_med,buffer),temp_ts,rep(temp_med,buffer))

    L = hampel_cust(temp_ts,k=buffer,t0=1)
  #  L = hampel_cust(temp_ts,k=buffer,t0=0.5)
    bool_arr = array(T,(act_len+(buffer*2)))
    bool_arr[L$ind] = F
    bool_arr = bool_arr[(buffer+1):(act_len+buffer)]
    red_bool[all_good] = bool_arr
  }

    ### extra second filter to remove obvious outliers from the red_vals set
    mn = mean(red_vals[red_bool])
    sd = sd(red_vals[red_bool],na.rm=T)
    thres = 2

    out_ind = !is.na(red_vals) & red_vals>mn-thres*sd & red_vals<mn+thres*sd
    bad_bool = red_bool
    red_bool[!out_ind] = F

    plot_flag = 0
    if(plot_flag ==1) {
      ## for making some example plots
      pdf("example_pixel.pdf")
      bad_pts = red_vals[!bad_bool & !red_bool & all_good]
      years = round(as.numeric(names(bad_pts))/1000,0)
      temp_good = red_vals[bad_bool & !red_bool]
      good_years = round(as.numeric(names(temp_good))/1000,0)
      plot(years,bad_pts,pch="*",col=2,cex=0.8,xlab="Time",ylab="Red Reflectance")
      points(good_years,temp_good,pch="*",col="orange",cex=0.8)
      temp_good2 = red_vals[red_bool]
      good_years2 = round(as.numeric(names(temp_good2))/1000,0)
      points(good_years2,temp_good2,pch="*",col=1,cex=0.8)
      text(paste("orig=",act_len,", final=",length(temp_good2),sep=""),x=1995,y=0.8)
      dev.off()
    }

  return(red_bool)
}

get_which_max <- function(x) {
  if(length(x[!is.na(x)])>0) {
    return(which(x == max(x, na.rm=T))[1])
  } else {
    return(NA)
  }
}


count_nas <- function(in_dat) {
  return(length(in_dat[!is.na(in_dat)]))
}

