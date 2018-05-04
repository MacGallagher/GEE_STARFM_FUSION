######## R code for GEE/STARFM ######### 
#### based on code by Faye Peters
#### modified by Megan Gallagher with help from Jake Graham
######################################################
### Used in conjunction with GEE code modified by Megan Gallagher, this is for Landsat8 and MODIS Terra Daily 250 m use
### Set up a folder that contains the Landsat, and MODIS files, as well as the starfm.exe and starfmconfig.txt
### create a separate folder called output in that folder
### for code changes and inquiries for adjustment for certain factors please email: megangallagher@u.boisestate.edu 

setwd("~/GEESTARFM")

library(raster)
library(rgdal)
library(zoo)
library(ggplot2)
library(fields)

###########################################Inputs##########################################################
## Load images and dates (be careful ofmatching dates and layers, the first layer may repeat)
## Read imagery dates and matching data:
image_dates<- read.csv("./2016_Dates_BOP", stringsAsFactors=F) 
## dates are not used in this code but are useful for other programs like bfastspatial, to check date overlap, and timing

modis<-brick("./2016_mod.tif")
landsat<-brick("./2016_landsat.tif")

##percentage of cover that is used to find "good" landsat scenes for interpolation 
## is the inverse of the actual cloud cover, 99 percent means 1 percent of the pixels is an actual value
perc<- 99

##################################################Code########################################################
#### STARFM code and function
starfm<-function(modis,landsat,perc){
config <- readLines("./StarFM_config.txt")
config <- gsub("(.*NROWS = ).*$", paste0("\\1", nrow(landsat)),config )
config <- gsub("(.*NCOLS = ).*$", paste0("\\1", ncol(landsat)),config )
cat(config, file="./StarFM_config.txt", sep="\n")

## Fix any missing data that enters as zeroes to our NA value.We have
## to do this by layer as operating on the entire stack may run into
## memory issues on larger subsets:
for (i in 1:nlayers(modis)) {
  ## Use the raster Which() function for speed:
  masked <- Which(modis[[i]] == 0, cells=TRUE)
  modis[[ i ]][ masked ] <- -32768
  masked <- Which(landsat[[i]] == 0, cells=TRUE)
  landsat[[ i ]][ masked ] <- -32768
}


## Automatically choose "good" landsat layers, at the moment this includes all actual landsat images
## if you want to create a threshhold for masking change the total percent of the pixel.
## i.e. if percent = 30, more than 70 percent of the image must be actual pixel values


area<-landsat@nrows*landsat@ncols
perc_Area<-(perc/100)*area

test2<-rep(NA,nlayers(modis))
filternew<-rep(NA,nlayers(modis))

for (i in 2:nlayers(modis)){
  test2[i] <-sum(landsat[[i]][]) 
  filternew[i] <-(test2[i]>{-32768*perc_Area}) == 1
}
filt3 = which(filternew==1)
filt3 <-append(filt3,1,0)

## Test landsat image from list

good_layer <-filt3(3)
plot(landsat[[good_layer]])

## If the above all works, then we run the following to loop over the
## MODIS time steps, filling in Landsat output as we go:

landsat_sim <- stack(modis)
landsat_sim[] <- NA

## Iterate and run StarFM for each MODIS date, choosing the
## nearest pair of good MODIS/Landsat dates, one before and
## one after the date being simulated where possible: 

good_landsat <- c(filt3) 
back <- 0 
forward <- 0
pb <- pbCreate((nlayers(landsat_sim)), "window", style=3,label='Time Step Progress')
for (i in 1:nlayers(landsat_sim)) {
  ##for (i in 1:10) {
  
  ## jakes/megan's work
  rm(back,forward)
  
  if (i %in% c(filt3)){
    back <- 0
    forward <- 0
    ls_t1temp<- i
    ls_t3temp <- i
  }
  else {
    foo <- good_landsat - i
    
    
    if(!length(which(foo==0))){
      back <- good_landsat[which(foo==max(foo[foo < 0]))]
      forward <- good_landsat[which(foo==min(foo[foo > 0]))]
      ls_t1temp <- back
      ls_t3temp <-forward
    }
  }
  
  ##end
  ls_t1 <- ls_t1temp
  ls_t3 <- ls_t3temp
  m_t1 <- ls_t1
  m_t3 <- ls_t3
  
  modis_t1 <- modis[[m_t1]]
  modis_t2 <- modis[[i]]
  modis_t3 <- modis[[m_t3]]
  landsat_t1 <- landsat[[ls_t1]]
  landsat_t3 <- landsat[[ls_t3]]
  
  
  writeRaster(modis_t1, filename="./modis_t1.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
  
  writeRaster(modis_t2, filename="./modis_t2.envi",bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
  
  writeRaster(modis_t3, filename="./modis_t3.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
  writeRaster(landsat_t1, filename="./landsat_t1.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
  writeRaster(landsat_t3, filename="./landsat_t3.envi", bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)
  
  system2(command="./StarFM.exe",args="StarFM_config.txt", wait=TRUE)
  
  landsat_t2_sim <- raster("./landsat_t2_sim.envi")
  
  
  ## Set any -32768 to NA values before writing:
  landsat_t2_sim[ landsat_t2_sim == -32768 ] <- NA
  landsat_sim[[i]] <- landsat_t2_sim[]
  
  ## In our filled data set, set any missing Landsat pixels to those
  ## simulated via StarFM:
  
  
    pbStep(pb, step=NULL, label='Processed Layer') }

   pbClose(pb, timer=T)
  return(landsat_sim)
}


landsat_sim<-starfm(modis,landsat,perc)


#### Possible raster outputs for saving
writeRaster(landsat_sim, filename="./output/2016_fusion.envi",
            bandorder='BSQ', datatype='INT2S', format="ENVI", overwrite=TRUE)

writeRaster(landsat_sim, filename="./output/2016_fusion.tif", bandorder='BSQ', 
            datatype='INT2S',format='GTiff', overwrite=TRUE)

#####################################################################################################################
#### Basic Gaussian Smoothing
######Inputs
upper<-8000 ##anyting above this counts as noise
lower<-1000 ## anything below this counts as noise
SD<-10 ## number of days for smoothing window
landsat_sim<-landsat_sim ## name of raster brick
############################################################
#gaussian kernel
GausKern <- function(SD){
  gk <- rep(NA, SD*7)
  inds <- seq(-(length(gk)/2), (length(gk)/2), 1) + .5
  for(i in 1:length(gk)){
    gk[i] <- (1/(sqrt(2*pi*SD^2))) * exp(-((inds[i])^2)/(2*SD^2))
  }
  gk <- gk + ((1-sum(gk))/length(gk))
  return(gk)
}
gk <- GausKern(SD)

# prepare smoothing function
# NAs are filled with median of that pixel through the time series (check bounds for best na)
eles <- landsat_sim@ncols * landsat_sim@nrows

for (i in 1:eles){
  
  tmp <- as.numeric(landsat_sim[i])
  tmp[tmp>upper] <-NA
  tmp[tmp<lower] <-NA
  tmp[is.na(tmp)] <- median(tmp,na.rm = TRUE)
  test[i] <- convolve(tmp[2:(length(tmp)-2)],gk, type = "filter")
}

writeRaster(test, filename="./output/landsat_smooth.tif", bandorder='BSQ', datatype='INT2S',format='GTiff', overwrite=TRUE)

###############################################Plots for smoothing view#############################################
#par(mfrow=c(2,1))
plot(1:length(landsat_sim[1]), landsat_sim[1], type = "l")
lines(1:length(landsat_sim[1]), test[1], col = "red", lwd = 2)



