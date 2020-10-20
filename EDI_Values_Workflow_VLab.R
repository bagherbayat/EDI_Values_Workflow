### EDI_Values_Workflow:
#This workflow quantifies one decade of (agricultural) water stress levels across Europe using satellite-derived Evapotraspiration (ET) data sets and Evaporative Drought Index values

## Authors: Bagher Bayat (b.bayat@fz-juelich.de and bagher.bayat@gmail.com) and Carsten Montzka (c.montzka@fz-juelich.de)
#Institute of Bio- and Geosciences: Agrosphere (IBG-3), Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
#Date:  18 March 2019

## Changelog:
#-system commands are used for reprojecting the input data sets (20 August 2020)
#-Implemented system commands are automatically executed for each file (20 August 2020)
#-The EDI outputs are better visalized (20 August 2020)
#-Spatial resolution (4 km) and resampling method (bilinear) are set in system commands (10 October 2020)
#-EDI workflow is running for a time-series from 2011-2020 (10 October 2020)

## Main inputs: 
#1. Time series of actual evapotranspiration (ETa) data set at daily step [mm] from 2011-2020 derived from the Spinning Enhanced Visible and Infrared Imager (SEVIRI) sensor onboard the Meteosat Second Generation (MSG) satellites 
#2. Time series of reference evapotranspiration (ET0) data set at daily step [mm] from 2011-2020 derived from the Spinning Enhanced Visible and Infrared Imager (SEVIRI) sensor onboard the Meteosat Second Generation (MSG) satellites
#3. Study area border (Europe) as a polygon shapefile

## Main outputs:
#1. Individual daily maps of water stress levels archived in a zip file
#2. Individual text reports (tables) containing water stress levels based on the percentage of the total land area archived in a zip file
#3. Time series of daily maps of water stress levels as an animated gif file
#4. Time series of daily reports (tables) containing water stress levels in a csv file

## Extent:
#European Union (can also be easily adapted for any specific country or Global)
#Spatial resolution: 4km
#Temporal resolution: daily (this can be adopted to weekly, monthly and yearly)

## Targeted Policy and indicator:
#SDG 6.4 (indicator 6.4.2: Levels of water stress)
#This workflow is developed within the European Commission HORIZON 2020 Program ERA-PLANET/GEOEssential project [grant number: 689443].

## Main reference:
#(Yao et al., 2010)


###################################################################################################
## 1. Load required packages----
library(sp)
library(raster)
library(prettymapr)
library(magick)
library(gifski)
library(xlsx)
library(zip)


## 2. Set Working directory
dir <- "./"            #Setting dir for VLab
setwd(dir)


## 4.2 Read  ETa data part 2 (to read from Disk region)
sdir <- "./ETa/" #Set working directory
list.filenames_ETa_Disk <-
  list.files(path = sdir, pattern = "Disk") # create a list from all HDF files of Disk region from the current directory
list.data_ETa_Disk <-
  list() #create an empty list that will serve as a container to receive the incoming files

for (i in 1:length(list.filenames_ETa_Disk))
  #create a loop to read all data sets
{
  print(paste(
    "Step 2: Processing Disk region ETa data set",
    i,
    "of",
    length(list.filenames_ETa_Disk)
  ))
  
  # Reprojecting the ETa data element (ET) of HDF files
  system(
    paste(
      'gdal_translate -a_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8 +no_defs" -a_ullr -5568000 5568000 5568000 -5568000 HDF5:',
      sdir,
      list.filenames_ETa_Disk[[i]],
      '://ET temp_DMET.tif',
      sep = ""
    )
  )
  
  system(
    paste(
      'gdalwarp -t_srs EPSG:4326 -te -10 33 34 73 -tr 0.04 0.04 -r bilinear -wo SOURCE_EXTRA=100 -overwrite temp_DMET.tif ET.tif',
      sep = ""
    )
  )
  
  
  # Read Reprojected file and apply the scalling
  setwd(dir)
  list.data_ETa_Disk[[i]] <- raster(paste(dir, "/ET.tif", sep = ""))
  list.data_ETa_Disk[[i]] <-
    list.data_ETa_Disk[[i]] / 1000           #scaling
  list.data_ETa_Disk[[i]][list.data_ETa_Disk[[i]] < 0] <-
    NA        #Excluding NEGATIVE values (In ET0 case it is -8000 to indicate missing values) in reference ET
  names(list.data_ETa_Disk[[i]]) <-
    list.filenames_ETa_Disk[[i]] #add the names of your data to the list
}


  



## 4. Read  ETa data
## 4.1 Read  ETa data part 1 (to read from Euro region)
sdir <- "./ETa/" #Set working directory
list.filenames_ETa_Euro <-
  list.files(path = sdir, pattern = "Euro") # create a list from all HDF files of Euro region from the current directory
list.data_ETa_Euro <-
  list() #create an empty list that will serve as a container to receive the incoming files

for (i in 1:length(list.filenames_ETa_Euro))
  #create a loop to read all data sets
{
  print(paste(
    "Step 2: Processing Euro region ETa data set",
    i,
    "of",
    length(list.filenames_ETa_Euro)
  ))
  
  # Reprojecting the ETa data element (ET) of HDF files
  system(
    paste(
      'gdal_translate -a_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8 +no_defs" -a_ullr -922623.5 5417891 4178899 3469966 HDF5:',
      sdir,
      list.filenames_ETa_Euro[[i]],
      '://ET temp_DMET.tif',
      sep = ""
    )
  )
  
  system(
    paste(
      'gdalwarp -t_srs EPSG:4326 -te -10 33 34 73 -tr 0.04 0.04 -r bilinear -wo SOURCE_EXTRA=100 -overwrite temp_DMET.tif ET.tif',
      sep = ""
    )
  )
  
  
  # Read Reprojected file and apply the scalling
  setwd(dir)
  list.data_ETa_Euro[[i]] <- raster(paste(dir, "/ET.tif", sep = ""))
  list.data_ETa_Euro[[i]] <-
    list.data_ETa_Euro[[i]] / 1000           #scaling
  list.data_ETa_Euro[[i]][list.data_ETa_Euro[[i]] < 0] <-
    NA        #Excluding NEGATIVE values (In ET0 case it is -8000 to indicate missing values) in reference ET
  names(list.data_ETa_Euro[[i]]) <-
    list.filenames_ETa_Euro[[i]] #add the names of your data to the list
}




# Combine two lists (list.data_ETa_Euro and list.data_ETa_Disk)
list.data_ETa <-
  do.call(c, list(list.data_ETa_Euro, list.data_ETa_Disk))

