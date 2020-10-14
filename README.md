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
