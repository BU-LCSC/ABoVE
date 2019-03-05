#######################################################
## Code - combines the output phenology/greenness metrics
## Takes 5 arguments:
##         Rscript make_tiffs.R in_tile in_dir out_dir cur_met cur_year
##
## in_tile - the ABOVE tile Bh##v##
## in_dir - the location of the chunked RData outputs from Phen/Green 
## out_dir - the location of the saved output TIFF files
##   cur_met - needs to be either tm5, etm, or phen - chooses the type of output
##  cur_year - only needed if the choices are tm5 or etm - gives the year out
## 
##
###  Last modified DSM - 1/24/18
##########################################################

## loads the packages we need mostly raster
require(raster)
library(RColorBrewer)
#library(doParallel)
#library(foreach)

## function to make a raster given a tile id starting with "B"
## cur_dat needs to be 6000x6000 pixels
make_ras <- function(cur_dat,tile) {
    ## ulx,uly of entire ABOVE grid
    ulx=-3400020
    uly=4640000
    ## pixel size and nrow/ncol
    pix=30
    dim=6000
    
    ## split up the tile id into x,y
    cur_tile_x = as.numeric(substr(tile,3,4))
    cur_tile_y = as.numeric(substr(tile,6,7))
    
    ## calculate the 4 corners of the map 
    ## given the upper left corner of the grid
    ulx_map= ulx + (cur_tile_x * pix * dim)
    uly_map= uly - (cur_tile_y * pix * dim)
    
    lrx_map= ulx + ((cur_tile_x+1) * pix * dim)
    lry_map= uly - ((cur_tile_y+1) * pix * dim)
    
    ## turn the data into a 6000x6000 matrix
    xy = matrix(as.numeric(cur_dat),nrow=dim,ncol=dim,byrow=T)
    # Turn the matrix into a raster
    rast <- raster(xy)
    # Give it x/y coords
    extent(rast) <- c(ulx_map,lrx_map,lry_map,uly_map)
    # ... and assign a projection
    projection(rast) <- CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
  
    return(rast)
}

## older function for making a JPG from the data
make_jpg <- function(ras,cur_name,out_file) {
    mycolors = colorRampPalette(c("red","orange","yellow","green","darkgreen","blue"))

    cur_range = round(quantile(in_dat,c(0.05,0.95),na.rm=T),2)
    cur_int = (cur_range[2]-cur_range[1])/5
    
      jpeg(file=out_file,width=5,height=5,units="in",res=250)
      par(mar=rep(0.2, 4), oma=rep(0.1, 4))
      #image(ras,col=mycolors(150))
    
    ## add a label with the year and date
    text(x=3000,y=500, labels=cur_name, col=1, cex=4)
    plot(ras, legend.only=TRUE, col=mycolors(150),
         legend.width = 2,
         axis.args=list(at=seq(cur_range[1], cur_range[2],cur_int),
                        labels=seq(cur_range[1],cur_range[2], cur_int), 
                        cex.axis=0.6),
         legend.args=list(text="          ", side=4, font=2, line=2.5, cex=0.8))
         #legend.args=list(text='NDVI', side=4, font=2, line=2.5, cex=0.8))
    dev.off()

}  ## end function


## function to write the pheno output to tiff files
## will be saved out as 68 tiff files
write_pheno_out <- function(in_tile,in_dir,out_dir) {

    fill = -1
    ## find all the files - each has a 2 row chunk of output
    in_files <- list.files(path=in_dir,pattern=glob2rx("*.RData"),full.names=T,include.dirs=T)
    num_files = length(in_files)

    ## setup the various parameters
    ## each chunk has 6000*2 pixels
    nrows=6000
    ncols=6000
    nrows_per_chunk = 2
    tot_pix_chunk = nrows_per_chunk*ncols

    ## make a list to hold all of the pheno metrics 68 x npix
    npix= nrows*ncols
    num_pheno = 68
    all_pheno_out = vector("list",num_pheno)
    for(b in 1:num_pheno) {
      all_pheno_out[[b]] = array(NA,npix)
    }

    time_process <- system.time  (
      
      for(i in seq(1,(nrows/nrows_per_chunk))) {
        #for(i in 1:100) {
        in_file = paste(in_dir,"/",in_tile,"_",i,".RData",sep="")
        if(file.exists(in_file)) {
          load(in_file)
        } else {
          print(paste("File ",in_file,"doesnt exist!"))
          next
        }
        start = ((i-1)*tot_pix_chunk)+1
        end = i*tot_pix_chunk
            ## load each chunk metric into a larger list item
           for(b in 1:num_pheno) {
                all_pheno_out[[b]][start:end] = pheno_out[,b]
        }
        
      }
    )  ## end system.time
    print(time_process)

    ## write out each to file 
    ## could write one file out or only a subset of the bands
    for(b in 1:num_pheno) {
      na_ind = is.na(all_pheno_out[[b]])
      all_pheno_out[[b]][na_ind] = fill
      pheno_ras = make_ras(all_pheno_out[[b]],in_tile)
      out_pheno_name = paste(out_dir,"/pheno.m",b,".tif",sep="")
      writeRaster(pheno_ras,out_pheno_name,format="GTiff",datatype='FLT4S',NAFlag=fill, overwrite=T)
    }
    
    ## doesnt return anything
}  ## end write pheno

## this function will output either tm5 or etm greenness metrics as tiffs
## the metrics it outputs are: red, nir, swir1, swir2, date of composite, number of obs
write_green_out <- function(in_tile,in_dir,out_dir,cur_met,y) {

    ## find all the files - each has a 2 row chunk of output
    in_files <- list.files(path=in_dir,pattern=glob2rx("*.RData"),full.names=T,include.dirs=T)
    num_files = length(in_files)

    ## setup the various parameters
    ## each chunk has 6000*2 pixels
    nrows=6000
    ncols=6000
    nrows_per_chunk = 2
    tot_pix_chunk = nrows_per_chunk*ncols

    npix= nrows*ncols

    ## only output the following bands - this would need to be changed in the composite_fxns.R code
    band_names = c("red","nir","swir1","swir2","num","date")
    num_bands = length(band_names)
    
    if(cur_met == "etm") {
      all_years = seq(1999,2014)
    } else {
      all_years = seq(1984,2011)
    }
    
    ### load values from arrays
    ## based on your year and desired output
    start_y = which(all_years==y)

    cur_out = array(NA,dim=c(npix,num_bands))

    time_process <- system.time  (

        for(i in seq(1,(nrows/nrows_per_chunk))) {
        #for(i in 1:100) {
          in_file = paste(in_dir,"/",in_tile,"_",i,".RData",sep="")
            if(file.exists(in_file)) {
              load(in_file)
            } else {
              print(paste("File ",in_file,"doesnt exist!"))
              next
            }
            ## start index and end index for the output array
            start = ((i-1)*tot_pix_chunk)+1
            end = i*tot_pix_chunk
            
            ## copy the data into the larger array for that year/chunk
            start_out = (start_y-1)*num_bands + 1
            end_out = start_y*num_bands
            ## inputs are saved in a large array with npix rows and nbands*nyears cols
            cur_out[start:end,] = green_out[[cur_met]][,start_out:end_out]

        }
        )  ## end system.time
      
        print(time_process)
      
      
      #cur_dat[cur_dat==fill]=NA
      ## make it a raster layer
      fill = -1
      na_ind = is.na(cur_out)
      ## fill missing with -1
      cur_out[na_ind] = fill
      ## make a stack of all the rasters
      cur_ras = make_ras(cur_out[,1],in_tile)
      for(b in 2:num_bands) {
        cur_ras = stack(cur_ras,make_ras(cur_out[,b],in_tile))
      }
      ## write the stack out to file
      out_name = paste(out_dir,"/",cur_met,".",y,".tif",sep="")
      writeRaster(cur_ras,out_name,format="GTiff",datatype='FLT4S',NAFlag=fill, overwrite=T)
} ## end write greenness


######################
#### BEGIN MAIN
###############

## read in arguments if debug flag is 0
debug=0
if(debug==0) {
  args = commandArgs(trailingOnly=T)
  in_tile = args[1]
  in_dir = args[2]
  out_dir = args[3]
  cur_met = args[4]
  cur_year = args[5]

} else {
  in_tile = "Bh04v06"
  in_dir = paste("../", in_tile,"_phen",sep="")
 
  out_dir=paste("./tifs_",in_tile,sep="")
  cur_met = "tm5"
  cur_year = 2000
}

## if the cur_met is not one of the 3 then abort
if(cur_met != "tm5" && cur_met != "etm" && cur_met != "phen") {
    print(paste("Wrong output selected",cur_met,"abort!"))
    quit(save="no")
}

## depending on the out metric we choose different paths here
if(cur_met=="phen") {
	print(paste("Writing out ",in_tile,in_dir,out_dir))
    write_pheno_out(in_tile,in_dir,out_dir)
## not pheno
} else {
    ## at this point it can only be tm5 or etm
    write_green_out(in_tile,in_dir,out_dir,cur_met,cur_year)
}

print("Done with all output!")





