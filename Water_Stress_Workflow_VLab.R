## Water Stress Workflow:
#  This is a simple script to quantify water stress by computing daily Evaporative Drought Index(EDI).
#  This index is also called the crop coefficient (Kc) in the literature.
#  In this script we use both of these terms interchangeably.

#  Authors: Bagher Bayat (b.bayat@fz-juelich.de/bagher.bayat@gmail.com), Carsten Montzka (c.montzka@fz-juelich.de)
#  Update:  18 March 2019

## Main inputs:
# 1. Actual evapotranspiration (ETa) at daily step [mm day-1]
# 2. Reference evapotranspiration (ET0) at daily step [mm day-1]
# 3. Projection files (lat and lon) [deg]
# 4. European (EU) countries border file [deg]

## Main outputs:
# Maps of water stress (drought) levels 
# Text reports (tables) containing various water stress levels for each country based on the percentage of the total land area 

## Spatial Extent:
# European Union (can also be easily adapted for any specific country or Global)

## Targeted Policy and indicator:
# SDG 6.4 (indicator 6.4.2: Levels of water stress)

#  Useful references:
#(Anderson et al., 2016; Bayat et al., 2018b, 2018a; Kim and Rhee, 2016; Narasimhan and Srinivasan, 2005)

###################################################################################################
# 1. load required packages----

require(sf)
library(sp)
library(rgdal)
library(raster)
library(h5)
library(Rcpp)
library(maptools)
library(rgeos)


# 2. Set Working directory

#dir<-"C:/Users/b.bayat/test_March/" #setting directory locally

dir<-"./"            #setting dir for VLab


# 3. Load global dataset

## 3.1. Read projection files (lat and lon)
setwd(dir)  
lat<-raster(paste(dir,"/projectionFiles/hdf5_lsasaf_msg_lat_msg-disk_4bytesprecision",sep=""))
lon<-raster(paste(dir,"/projectionFiles/hdf5_lsasaf_msg_lon_msg-disk_4bytesprecision",sep=""))

lat <- lat[]/10000  #Scaling data
lon <- lon[]/10000  #Scaling data
valid <- (lon < 35  & lon > -10) & (lat < 75 & lat > 35) #Extracting Europe area from the files

lat <- lat[valid]
lon <- lon[valid]


## 3.2. Read ET0 data

sdir<-"./ET0/" #set working directory
list.filenames_ET0<-list.files(path=sdir,pattern="HDF") #list all HDF files from the current directory
list.data_ET0<-list() #create an empty list that will serve as a container to receive the incoming files


for (i in 1:length(list.filenames_ET0))  #create a loop to read in your data
{
  
  list.data_ET0[[i]]<-h5file(paste(sdir,list.filenames_ET0[i],sep=""))    # reads HDF files
  list.data_ET0[[i]]<-list.data_ET0[[i]]["/METREF"]    # opens "/METREF" element
  list.data_ET0[[i]]<-readDataSet(list.data_ET0[[i]])  # reads "/METREF" values as matrix
  list.data_ET0[[i]]<-list.data_ET0[[i]]/100           #scaling
  list.data_ET0[[i]]<-raster(list.data_ET0[[i]])       #converting to a raster
  list.data_ET0[[i]] <- list.data_ET0[[i]][valid]      # croping European area
  list.data_ET0[[i]][list.data_ET0[[i]]<0] <- NA       #Excluding NEGATIVE values (gaps) in reference ET
}

names(list.data_ET0)<-list.filenames_ET0 #add the names of your data to the list


## 3.3. Read ETa data

sdir<-"./ETa/"  #set working directory
list.filenames_ETa<-list.files(path=paste(sdir),pattern="HDF") #list all HDF files from the current directory
list.data_ETa<-list() #create an empty list that will serve as a container to receive the incoming files


for (j in 1:length(list.filenames_ETa)) #create a loop to read in your data
  
{
  list.data_ETa[[j]]<-h5file(paste(sdir,list.filenames_ETa[i],sep=""))    # reads HDF files
  list.data_ETa[[j]]<-list.data_ETa[[j]]["/ET"]    # opens "/ET" element
  list.data_ETa[[j]]<-readDataSet(list.data_ETa[[j]])  # reads "/ET" values as matrix
  list.data_ETa[[j]]<-list.data_ETa[[j]]/1000           #scaling
  list.data_ETa[[j]]<-raster(list.data_ETa[[j]])       #converting to a raster
  list.data_ETa[[j]] <- list.data_ETa[[j]][valid]      # croping Europe area
  list.data_ETa[[j]][list.data_ETa[[j]]<0] <- NA  #Excluding NEGATIVE values (gaps) in reference ET
}

names(list.data_ETa)<-list.filenames_ETa #add the names of your data to the list


# 4. Compute crop coefficient (Kc) (This is also called Evaporative Drought Index (EDI))

list.data_Kc<-list()# 3.1. create an empty list that will serve as a container to receive the incoming files

## 4.1. create a loop to read in your data

for(i in 1:length(list.filenames_ETa)){
  list.data_Kc[[i]]<- 1-(list.data_ETa[[i]]/ list.data_ET0[[i]])   # computing Kc
  
  ## 4.2. Adding projection
  #### 4.2.1. defining an SpatialPointsDataFrame object using the lat and lon as the coordinates, the ET file as the data, and with the long/lat projection.
  
  projLL <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
  list.data_Kc[[i]] <- SpatialPointsDataFrame(cbind(lon, lat), data=data.frame(list.data_Kc[[i]]),proj4string=projLL)
  
  #### 4.2.2. Reprojectting this SpatialPointsDataFrame to the GEOS projection with spTransform
  
  projGeos <- CRS("+proj=geos +lon_0=0 +h=35785831 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs ")
  xy <- spTransform(list.data_Kc[[i]], projGeos)
  
  #### 4.2.3.Defining an SpatialPixelsDataFrame (a regular grid) with the new coordinates, which can be easily converted to a RasterLayer object.
  
  list.data_Kc[[i]] <- SpatialPixelsDataFrame(xy, data.frame(list.data_Kc[[i]]), tolerance=0.311898)
  list.data_Kc[[i]] <- raster(list.data_Kc[[i]])
  
  
  # 4.3. Set LB and UB for Kc 
  list.data_Kc[[i]][list.data_Kc[[i]]<0] <- 0    #Minimum Kc is 0 (LB)
  list.data_Kc[[i]][list.data_Kc[[i]]>1] <- 1    #Maximum Kc is 1 (HB)
  
  
  # 4.4. Saving and Plotting daily time series of Kc
  
  ## 4.4.1. trying to get names from input files
  
  list.filenames_ETa[i]<-strsplit(list.filenames_ETa[i],split='Disk_')[[1]][2]  #spliting the title to two parts and taking the date
  list.filenames_ETa[i]<- gsub('.{4}$', '', list.filenames_ETa[i])
  list.filenames_ETa[i]<- paste("Water Stress Levels_", list.filenames_ETa[i], sep=" ")
  
  
  ## 4.4.2. masking the map based on European countries border
  
  ### This runs locally 
  ##sdir<-"./EU_Border/" #set working directory
  ##unzip(zipfile="./EU_Border/data.zip",exdir="./EU_Border/data")#unzipping the data folder
  ##file<-paste(sdir,"/data/NUTS_RG_01M_2013.shp",sep="")  
  ##europe.map<- shapefile(file) #reading unzipped shapefile
  
  ### This runs on VLab 
  sdir<-"./EU_Border/" #set working directory
  system("unzip ./SHP/data.zip -d ./SHP/")
  file<-paste(dir,"/SHP/NUTS_RG_01M_2013.shp",sep="")  
  europe.map<- shapefile(file) 
  
  europe.map <- europe.map[europe.map$STAT_LEVL_ == 0,] #reading country (state) level data 
  
  project2<-"+proj=geos +lon_0=0 +h=35785831 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" #adding projection to the shapefile
  europe.map <-spTransform(europe.map, project2)
  
  list.data_Kc[[i]] <- mask(x = list.data_Kc[[i]], mask = europe.map) #masking Kc (EDI) products
  
  
  # 4.5. classifying the Kc map
  
  ### 4.5.1. create classification matrix
  reclass_df <- c(0, 0.225, 0,
                  0.225, 0.475, 1,
                  0.475, 0.725, 2,
                  0.725, 0.975, 3,
                  0.975, Inf, 4)
  
  ### 4.5.2. reshape the object into a matrix with columns and rows
  reclass_m <- matrix(reclass_df,
                      ncol = 3,
                      byrow = TRUE)
  
  ### 4.5.3. reclassify the raster using the reclass object - reclass_m
  list.data_Kc[[i]] <- reclassify(list.data_Kc[[i]],reclass_m)
  
  ### 4.5.4. plot reclassified data
  
  r_colors <- c("yellow", "burlywood1", "darkgoldenrod1","red", "darkred")
  plot(list.data_Kc[[i]],
       legend = FALSE,
       col = r_colors)
  
  legend("topright",
         legend = c("Near Normal","Moderate Stress", "Severe Stress", "Extreme Stress", "Exceptional Stress"),
         fill = r_colors,
         border = FALSE,
         bty = "n") # turn off legend border
  plot(europe.map, add=TRUE,lwd = 0.2)
  title(list.filenames_ETa[i])
  
  
  # 5. Saving daily time series as jpg files
  
  #jpeg(paste(dir,list.filenames_ETa[i],'.jpg',sep=""),width=725, height=512)  # the file name is selected based on input file date
  jpeg(paste(dir,"Water_Stress_Map.jpg",sep=""),width=725, height=512)                # the file name is fixed as "Water_Stress_Map" (underlines [_] are important for VLab)
  
  
  plot(list.data_Kc[[i]],
       legend = FALSE,
       col = r_colors)
  
  legend("topright",
         legend = c("Near Normal","Moderate Stress", "Severe Stress", "Extreme Stress", "Exceptional Stress"),
         fill = r_colors,
         border = FALSE,
         bty = "n") # turn off legend border
  plot(europe.map, add=TRUE,lwd = 0.1)
  title(list.filenames_ETa[i])
  dev.off()
  
  # 6. Saving daily time series as txt reports
  
  ## 6.1. Extract raster values to polygons                             
  (v <- extract(list.data_Kc[[i]], europe.map[1]))
  
  
  ## 6.2. Get class counts for each polygon
  v.counts <- lapply(v,table)
  
  ## 6.3. Calculate class percentages for each polygon
  (v.pct <- lapply(v.counts, FUN=function(x){ (x / sum(x))*100 }))
  
  ## 6.4. Seperate columns of v.pct
  out0 <- lapply(  v.pct , '[', '0')
  out1 <- lapply(  v.pct , '[', '1')
  out2 <- lapply(  v.pct , '[', '2')
  out3 <- lapply(  v.pct , '[', '3')
  out4 <- lapply(  v.pct , '[', '4')
  
  ## 6.5.Replace missing values with NA in different columns
  l0 <- sapply( out0 , function(x) c( x , rep( NA , 1-length(x) ) ) ) 
  l1 <- sapply( out1 , function(x) c( x , rep( NA , 1-length(x) ) ) ) 
  l2 <- sapply( out2 , function(x) c( x , rep( NA , 1-length(x) ) ) ) 
  l3 <- sapply( out3 , function(x) c( x , rep( NA , 1-length(x) ) ) ) 
  l4 <- sapply( out4 , function(x) c( x , rep( NA , 1-length(x) ) ) ) 
  
  
  NearNormal <-lapply(l0, unlist)
  ModerateStress<-lapply(l1, unlist)
  SevereStress<-lapply(l2, unlist)
  ExtremeStress<-lapply(l3, unlist)
  ExceptionalStress<-lapply(l4, unlist)
  
  
  ## 6.6. creating an empty datafram
  df<-as.data.frame(matrix(ncol=5,nrow=35))
  names(df)<-c(1,2,3,4,5)
  
  ## 6.7. Fill the empty datafram with columns
  
  df[1]<-do.call(rbind, NearNormal)
  df[2]<-do.call(rbind, ModerateStress)
  df[3]<-do.call(rbind, SevereStress)
  df[4]<-do.call(rbind, ExtremeStress)
  df[5]<-do.call(rbind, ExceptionalStress)
  
  class.df<-df
  
  
  ## 6.8. Replace NA's with 0 and add names and naming the columns
  class.df[is.na(class.df)] <- 0 
  colnames(class.df) <- c("NearNormal","ModerateStress","SevereStress","ExtremeStress", "ExceptionalStress")
  
  
  ## 6.9. adding a new columns for EU countries names
  
  class.df$EUcountry=europe.map[[1]]
  class.df[c("EUcountry","NearNormal","ModerateStress","SevereStress","ExtremeStress", "ExceptionalStress")]
  
  ## 6.7. Here is the txt report
  
  #txtfile <-  paste(dir,list.filenames_ETa[i],'.txt',sep="") # the file name is selected based on input file date
  txtfile <-  paste(dir,"Water_Stress_Report.txt",sep="") # the file name is fixed as "Water_Stress_Report" (underlines [_] are important for VLab)
  write.table(class.df[c("EUcountry","NearNormal","ModerateStress","SevereStress","ExtremeStress", "ExceptionalStress")], file =  txtfile, sep = " ", row.names = F, col.names = T,
              quote = F)
  
}
