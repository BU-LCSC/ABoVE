## created to read in a chunk of cached data Landsat time series data in NetCDF

## to load composite_fxns.R
require(ncdf4)
# require(data.table)
library("foreach")
library("iterators")
library("doParallel")
library("raster")
#library("outliers")

## location of our homemade functions
in_lib="/discover/nobackup/projects/landsatts/dsm_work/all_dsm_codes/dsm/above/read_nc/"
source(paste(in_lib,"composite_fxns_v1.2.R",sep=""))
source(paste(in_lib,"pheno_fxns_v1.2.R",sep=""))

debug=0
if(debug==0) {
  args = commandArgs(trailingOnly=T)
  in_dir = args[1]
  in_tile = args[2]
  out_dir = args[3]
  dist_dir = args[4]
  start_row = as.numeric(args[5])
  end_row = as.numeric(args[6])

} else {
  in_dir = "/projectnb/landsat/users/dsm/above/ccdc_develop/ccdc_Bh11v11_from_nc/Bh11v11_NC"
  in_tile = "Bh11v11"
  out_dir="./Bh11v11_green_110217"
  dist_dir="../discover/alaska_mos_110217/alaska_dist"  
  start_row = 1
  end_row = 10
}

#in_files <- list.files(path=in_dir,pattern=glob2rx("*.nc"),full.names=F,include.dirs=F)

#num_files = length(in_files)

## do some fiddling with the file names - to get the listing the way we want it
## need to think of a better way to do this
nrow = 2
total_nrows=6000
ncols = total_nrows
ncores=8
nchunks = total_nrows/nrow

num_files = nchunks

dist_file = paste(dist_dir,"/",in_tile,".dates.tif",sep="")
  if(file.exists(dist_file)) {
         all_dist_dat = as.vector(raster(dist_file))
  } else {
        print(paste("No dist file found. Expecting:",dist_file))
  }

#Register the parallel backend
registerDoParallel(ncores)

## setup to loop through files
for(i in start_row:end_row) {
#for(i in 1:10) {
  start = (i-1)*nrow
  end = (i*nrow) - 1

  ## same number of columns and rows 
   ## this gives the starting and ending pixels
  start_pix = (start * ncols) + 1
  end_pix = ((end+1) * ncols)
  cur_dist_dat = all_dist_dat[start_pix:end_pix]
    

  cur_file = paste(in_dir,"/",start,"-",end,".nc",sep="")
  
    print(paste("Now processing row ",i," and file ", cur_file,sep=""))

    ## name of the output file
    out_name = paste(out_dir,"/",in_tile,"_",i,".RData",sep="")
    ## initialize the outputs in case there is an error
    green_out = NULL
    pheno_out = NULL

     ## only create output if it doesnt already exist - else skip
     ##if(file.exists(out_name)) {
      ##print(paste("Output already exists at this location ",out_name,sep=""))
      ##next
     ##}
  
      ## check if input exists - else skip
      if(file.exists(cur_file)) {
        in_dat = load_nc_file(cur_file)
      } else {
        print(paste("Input does not exist at this location ",cur_file,sep=""))
        next
      }
  
      ## now it is ready to go - each item in the out list is a band (1,2,3,4,5,7,6,fmask - or whichever you want)
      ## each col in the output is a single date and each row is a pixel
      ## the cols should correspond to the dates in head_dat - sorted
      ## to see the time series for a pixel - reference it as out_list[[band]][pix,]
      
      ref_list = in_dat[[1]]
      head_dat = in_dat[[2]]
      rm(in_dat)
      names(ref_list) = c("blue","red","nir","swir1","swir2")

      ## we calculate NDVI from the ref_list object
      ## only extracted blue (1), red (2), and nir (3) bands
      evi_dat <- 2.5*(ref_list[["nir"]]-ref_list[["red"]])/
              (ref_list[["nir"]]+(6*ref_list[["red"]])-(7.5*ref_list[["blue"]])+1)
      evi_dat[evi_dat>1 | evi_dat<0] = NA
      
      ## Extra cloud screening in case Fmask is poor
      red_dat <- ref_list[["red"]]
      evi_dat[red_dat>0.1] = NA      
  
      #In case we're using a large enough matrix, switch to data table format
      #evi_dat <- data.table(evi_dat,keep.rownames=FALSE)
      
      ### green out will contain the max ndvi values for each year in the series. Needs the red and NDVI matrices as inputs as well as the header dat.
      time_process = system.time( green_out <- comp_ndvi(ref_list,head_dat) )
      
      print("Time to composite ndvi:")
      print(time_process)

      ### pheno_out will contain the phenology observations for each year in the series.  Needs the EVI matrix.
      time_process = system.time( pheno_out <- run_pheno(evi_dat,head_dat,cur_dist_dat) )
      
      print("Time to calc pheno:")
      print(time_process)
      
      ## one alternative is to save the intermediate files - they are relatively small
      #save(evi_dat,ndvi_dat,head_dat,file=temp_name)
      save(pheno_out,green_out,head_dat,file=out_name)

} ## end for loop


