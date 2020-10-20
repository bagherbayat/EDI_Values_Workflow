### EDI_Values_Workflow:
#This workflow quantifies one decade of (agricultural) water stress levels across Europe using satellite-derived Evapotraspiration (ET) data sets and Evaporative Drought Index values

##Authors: Bagher Bayat (b.bayat@fz-juelich.de and bagher.bayat@gmail.com) and Carsten Montzka (c.montzka@fz-juelich.de)
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

## 3. Read  ET0 data
sdir <- "./ET0/" #set working directory
list.filenames_ET0 <-
  list.files(path = sdir, pattern = "Disk") # create a list from all HDF files from the current directory
list.data_ET0 <-
  list() #create an empty list that will serve as a container to receive the incoming files

for (i in 1:length(list.filenames_ET0))
  #create a loop to read all data sets
{
  print(paste(
    "Step 1: Processing ET0 data set",
    i,
    "of",
    length(list.filenames_ET0)
  ))
  
  # Reprojecting the ET0 data element (METREF) of HDF files
  system(
    paste(
      'gdal_translate -a_srs "+proj=geos +h=35785831 +a=6378169 +b=6356583.8 +no_defs" -a_ullr -5568000 5568000 5568000 -5568000 HDF5:',
      sdir,
      list.filenames_ET0[[i]],
      '://METREF temp_METREF.tif',
      sep = ""
    )
  )
  
  system(
    paste(
      'gdalwarp -t_srs EPSG:4326 -te -10 33 34 73 -tr 0.04 0.04 -r bilinear -wo SOURCE_EXTRA=100 -overwrite temp_METREF.tif METREF.tif',
      sep = ""
    )
  )
  
  
  # Read Reprojected file and apply the scalling
  setwd(dir)
  list.data_ET0[[i]] <- raster(paste(dir, "/METREF.tif", sep = ""))
  list.data_ET0[[i]] <- list.data_ET0[[i]] / 100           #scaling
  list.data_ET0[[i]][list.data_ET0[[i]] < 0] <-
    NA        #Excluding NEGATIVE values (In ET0 case it is -8000 to indicate missing values) in reference ET
  names(list.data_ET0[[i]]) <-
    list.filenames_ET0[[i]] #add the names of your data to the list
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


## 4.2 Read  ETa data part 2 (to read from Disk region)
sdir <- "./ETa/" #set working directory
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


# Combine two lists (list.data_ETa_Euro and list.data_ETa_Disk)
list.data_ETa <-
  do.call(c, list(list.data_ETa_Euro, list.data_ETa_Disk))

## 5. Compute Evaporative Drought Index (EDI)

list.data_EDI <-
  list()# create an empty list that will serve as a container to receive the incoming files

#We first there is a need to remove xlsx file in the working directory generated from previous runs (Overwriting an existing xlsx worksheet is problematic!)
xlsx_file <- list.files(path = dir, pattern = "xlsx")
unlink(paste(dir, xlsx_file, sep = ""))

#Lets continue with the loop
for (i in 1:length(list.filenames_ET0))
{
  print(paste(
    "Step 3: Computing EDI values",
    i,
    "of",
    length(list.filenames_ET0)
  ))
  
  list.data_EDI[[i]] <-
    1 - (list.data_ETa[[i]] / list.data_ET0[[i]])   # computing EDI
  
  # Set LB for EDI
  list.data_EDI[[i]][list.data_EDI[[i]] < 0] <-
    0    #Minimum EDI is 0 (LB)
  
  
  # Trying to get names from input files. For this it is better to use ET0 since these data sets have always fixed file naming, while, that is changing for ETa
  list.filenames_ET0[i] <-
    strsplit(list.filenames_ET0[i], split = 'Disk_')[[1]][2]  #spliting the title to two parts and taking the date
  list.filenames_ET0[i] <- gsub('.{4}$', '', list.filenames_ET0[i])
  list.filenames_ET0[i] <-
    paste("Water_Stress_Levels_", list.filenames_ET0[i], sep = " ")#(underlines [_] in file names are important for VLab)
  
  # 6. Masking the map based on European countries border
  
  ### This runs locally
  sdir <- "./EU_Border/" #set working directory
  unzip(zipfile = "./EU_Border/data.zip", exdir = "./EU_Border/data")#unzipping the data folder
  file <- paste(sdir, "/data/NUTS_RG_01M_2013_Update.shp", sep = "")
  europe.map <- shapefile(file) #reading unzipped shapefile
  
  # ### This runs on VLab
  # sdir<-"./EU_Border/" #set working directory
  # system("unzip ./SHP/data.zip -d ./SHP/")
  # file<-paste(dir,"/SHP/NUTS_RG_01M_2013_Update.shp",sep="")
  # europe.map<- shapefile(file)
  
  europe.map <-
    europe.map[europe.map$STAT_LEVL_ == 0,] #reading country (state) level data
  
  project <-
    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  europe.map <- spTransform(europe.map, project)
  e <- extent(-10, 34, 33, 73) #This is EU border extent
  europe.map <- crop(europe.map, e)
  
  list.data_EDI[[i]] <-
    mask(x = list.data_EDI[[i]], mask = europe.map) #masking EDI products
  
  
  ## 7. Classifying the EDI map
  #Create classification matrix
  reclass_df <-  c(-Inf, 0.2, 1,
                   0.2, 0.4, 2,
                   0.4, 0.6, 3,
                   0.6, 0.8, 4,
                   0.8, 1, 5,
                   1, Inf, 6)
  
  
  #Reshape the object into a matrix with columns and rows
  reclass_m <- matrix(reclass_df,
                      ncol = 3,
                      byrow = TRUE)
  
  #Reclassify the raster using the reclass object - reclass_m
  list.data_EDI[[i]] <-
    reclassify(list.data_EDI[[i]], reclass_m)
  
  #Plot reclassified data
  r_colors <-
    c("greenyellow",
      "yellow",
      "burlywood1",
      "darkgoldenrod1",
      "red",
      "darkred")
  
  #Margins for our plot
  par(mar = c(4, 4, 1.2, 0)) # Set the margin on all sides  # the syntax par(mar = c(bottom, left, top, right))
  
  plot(
    list.data_EDI[[i]],
    breaks = 0:6,
    xlim = c(-20, 40),
    ylim = c(30, 75),
    legend = FALSE,
    col = r_colors,
    xlab = "Longitude [deg]",
    ylab = "Latitude [deg]"
  ) #breaks are needed for assigning colors to classes correctly
  
  plot(
    europe.map,
    add = TRUE,
    lwd = 1,
    xlim = c(-20, 40),
    ylim = c(30, 75)
  ) # add shape to map
  
  addnortharrow(
    pos = "topright",
    padin = c(0.15, 0.15),
    scale = 0.65,
    lwd = 0.65,
    border = "black",
    cols = c("white", "black"),
    text.col = "black"
  )
  
  
  # Using raster library for scalebar
  raster::scalebar(
    d = 1000,
    # distance in km
    xy = c(-40, 34),
    type = "bar",
    divs = 4,
    below = "km",
    lonlat = TRUE,
    label = c(0, 500, 1000),
    adj = c(0, -0.75),
    lwd = 1
  )
  
  mtext(
    "GCS WGS 1984",
    side = 1,
    line = -1.43,
    cex = 0.9,
    at = 50
  )
  
  
  legend(
    "topleft",
    legend = c("<0.2",
               "0.2 - 0.4",
               "0.4 - 0.6",
               "0.6 - 0.8",
               "0.8 - 1",
               "> 1"),
    
    fill = r_colors,
    border = "black",
    bty = "n",
    # turn off legend border
    title = "EDI values [-]"
  )
  
  
  title(list.filenames_ET0[i])
  
  
  ## 8. Saving daily time series as jpg files
  
  dpi <- 500
  jpeg(
    paste(dir, list.filenames_ET0[i], '.jpg', sep = ""),
    width = 10 * dpi,
    height = 6 * dpi,
    res = dpi
  )  # the file name is selected based on input file date
  
  #Margins for our plot
  par(mar = c(4, 4, 1.2, 0)) # Set the margin on all sides  # the syntax par(mar = c(bottom, left, top, right))
  
  plot(
    list.data_EDI[[i]],
    breaks = 0:6,
    xlim = c(-20, 40),
    ylim = c(30, 75),
    legend = FALSE,
    col = r_colors,
    xlab = "Longitude [deg]",
    ylab = "Latitude [deg]"
  ) #breaks are needed for assigning colors to classes correctly
  
  plot(
    europe.map,
    add = TRUE,
    lwd = 1,
    xlim = c(-20, 40),
    ylim = c(30, 75)
  ) # add shape to map
  
  addnortharrow(
    pos = "topright",
    padin = c(0.15, 0.15),
    scale = 0.65,
    lwd = 0.65,
    border = "black",
    cols = c("white", "black"),
    text.col = "black"
  )
  
  
  # Using raster library for scalebar
  raster::scalebar(
    d = 1000,
    # distance in km
    xy = c(-40, 34),
    type = "bar",
    divs = 4,
    below = "km",
    lonlat = TRUE,
    label = c(0, 500, 1000),
    adj = c(0, -0.75),
    lwd = 1
  )
  
  mtext(
    "GCS WGS 1984",
    side = 1,
    line = -1.43,
    cex = 0.9,
    at = 50
  )
  
  
  legend(
    "topleft",
    legend = c("<0.2",
               "0.2 - 0.4",
               "0.4 - 0.6",
               "0.6 - 0.8",
               "0.8 - 1",
               "> 1"),
    
    fill = r_colors,
    border = "black",
    bty = "n",
    # turn off legend border
    title = "EDI values [-]"
  )
  
  
  title(list.filenames_ET0[i])
  dev.off()
  
  
  ## 9. Lets generate some statistics per country as a text report
  
  #Extract raster values to polygons
  (v <- extract(list.data_EDI[[i]], europe.map[1]))
  
  
  #Get class counts for each polygon
  v.counts <- lapply(v, table)
  
  #Calculate class percentages for each polygon
  (v.pct <- lapply(
    v.counts,
    FUN = function(x) {
      (x / sum(x)) * 100
    }
  ))
  
  #Seperate columns of v.pct
  out1 <- lapply(v.pct , '[', '1')
  out2 <- lapply(v.pct , '[', '2')
  out3 <- lapply(v.pct , '[', '3')
  out4 <- lapply(v.pct , '[', '4')
  out5 <- lapply(v.pct , '[', '5')
  out6 <- lapply(v.pct , '[', '6')
  
  
  #Replace missing values with NA in different columns
  l1 <-
    sapply(out1 , function(x)
      c(x , rep(NA , 1 - length(x))))
  l2 <-
    sapply(out2 , function(x)
      c(x , rep(NA , 1 - length(x))))
  l3 <-
    sapply(out3 , function(x)
      c(x , rep(NA , 1 - length(x))))
  l4 <-
    sapply(out4 , function(x)
      c(x , rep(NA , 1 - length(x))))
  l5 <-
    sapply(out5 , function(x)
      c(x , rep(NA , 1 - length(x))))
  l6 <-
    sapply(out6 , function(x)
      c(x , rep(NA , 1 - length(x))))
  
  
  Class_1 <- lapply(l1, unlist)
  Class_2 <- lapply(l2, unlist)
  Class_3 <- lapply(l3, unlist)
  Class_4 <- lapply(l4, unlist)
  Class_5 <- lapply(l5, unlist)
  Class_6 <- lapply(l6, unlist)
  
  
  #Creating an empty datafram
  df <- as.data.frame(matrix(ncol = 6, nrow = 34))
  names(df) <- c(1, 2, 3, 4, 5, 6)
  
  #Filling the empty datafram with columns
  df[1] <- do.call(rbind, Class_1)
  df[2] <- do.call(rbind, Class_2)
  df[3] <- do.call(rbind, Class_3)
  df[4] <- do.call(rbind, Class_4)
  df[5] <- do.call(rbind, Class_5)
  df[6] <- do.call(rbind, Class_6)
  
  class.df <- df
  
  
  #Replace NA's with 0 and add names and naming the columns
  class.df[is.na(class.df)] <-
    0    #This changes NA values to zero in the report for certain classes that do not have any percentage (for instance, CZ do not have any near normal, moderate classes, so they are inserted as NA in the report so we need to change this to zero
  colnames(class.df) <-
    c("[<0.2]",
      "[0.2 - 0.4]",
      "[0.4 - 0.6]",
      "[0.6 - 0.8]",
      "[0.8 - 1]",
      "[>1]")
  
  
  #Adding a new columns for EU countries names
  class.df$Country = europe.map[[1]]
  class.df[c("Country",
             "[<0.2]",
             "[0.2 - 0.4]",
             "[0.4 - 0.6]",
             "[0.6 - 0.8]",
             "[0.8 - 1]",
             "[>1]")]
  

  #Here is the individual daily csv reports
    csvfile <-
    paste(dir, list.filenames_ET0[i], '.csv', sep = "") # the file name is selected based on input file date
  write.table(
    class.df[c("Country",
               "[<0.2]",
               "[0.2 - 0.4]",
               "[0.4 - 0.6]",
               "[0.6 - 0.8]",
               "[0.8 - 1]",
               "[>1]")],
    file =  csvfile,
    sep = ",",
    row.names = F,
    col.names = T,
    quote = F
  )

  # #Here is all reports in one xls file
  xlsfile <-
  paste(dir, "TimeSeries_Reports.xlsx", '.xlsx', sep = "") 
    write.xlsx(
    class.df[c("Country",
               "[<0.2]",
               "[0.2 - 0.4]",
               "[0.4 - 0.6]",
               "[0.6 - 0.8]",
               "[0.8 - 1]",
               "[>1]")],
    file =  xlsfile,
    sheetName = list.filenames_ET0[i],
    append = TRUE
  )
 
}

#10. Making gif (movie) from time series of EDIs
setwd(dir)
files_jpg <- list.files(path = dir, pattern = "jpg")
frames <- image_read(paste(dir, files_jpg, sep = ""))
image_write_gif(frames, path = "TimeSeries_Maps.gif", delay = 1) #delay is the duration of each frame in seconds

#11. Making zip files to save individual (daily) outputs
zipfile_jpg <- "Individual_Maps.zip"
zip(zipfile_jpg,paste(dir, files_jpg, sep = ""),recurse = F,compression_level = 9,include_directories = F,root = ".", mode = c("cherry-pick"))

zipfile_csv <- "Individual_Reports.zip"
files_csv <- list.files(path = dir, pattern = "csv")
zip(zipfile_csv,paste(dir, files_csv, sep = ""),recurse = F,compression_level = 9,include_directories = F,root = ".", mode = c("cherry-pick"))


#12. Finally removing intermadiate files in the directory 
setwd(dir)
files <- list.files(path = dir, pattern = "jpg")
unlink(paste(dir, files, sep = ""))
files <- list.files(path = dir, pattern = "csv")
unlink(paste(dir, files, sep = ""))
files <- list.files(path = dir, pattern = "tif")
unlink(paste(dir, files, sep = ""))

























