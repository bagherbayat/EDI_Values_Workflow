# WaterStress

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
