#Download Rangeland Analysis Platform Vegetation 
#Date Modified: 3/17/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script downloads RAP data from Earth engine for Mesonet stations. 

#load required libraries
library(reticulate) # allows for python interfacing
library(rgee) # R wrapper for the python GEE library
library(sf) # simple feature library - used for vectors
library(tidyverse) # package for tidy syntax etc
library(geojsonio) # package to send ROI SF objects to GEE

#set up the gee environment
use_condaenv("gee-base", conda = "auto",required = TRUE)
ee = import("ee")
ee_Initialize(drive = TRUE)

stations = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/All_Data') %>%
  substr(., 17,24)

#import mesonet data and convert to SF object (replace here with shp file)

mesonet = read.csv('https://mesonet.climate.umt.edu/api/v2/stations/?type=csv&clean=true') %>%
  #filter(station %in% stations) %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs(4326))


#Import RAP Data
RAP = ee$ImageCollection("projects/rangeland-analysis-platform/vegetation-cover-v3")$
  filter(ee$Filter$date('2017-01-01', '2022-01-01'))$
  #filter(ee$Filter$dayOfYear(152, 273))$
  toBands()

#extract by points
point_extract = ee_extract(
 #RAP image
  x = RAP,
  #sf object (no need to convert to ee object now)
  y = mesonet,
  #sampling function - mean works, you are extracting at native res
  fun = ee$Reducer$max(),
  #define resolution for extraction (30m = native)
  scale = 30,
  #tell rgee the shp file is a sf object
  sf = T
) 


#RAP Extract 
RAP_data = point_extract %>%
  st_drop_geometry() %>%
  select(-c(date_installed, sub_network, elevation, mesowest_id, gwic_id)) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(station, name), names_to = 'date') %>%
  mutate(year = substr(date, 2, 5),
         land_cover = substr(date,7,9)) %>%
  select(station, name, value, year, land_cover)


#Calculate averages 
average_RAP = RAP_data %>%
  group_by(station, name, land_cover) %>%
  summarise(mean_percent = mean(value))

RAP_final = as.data.frame(average_RAP)


write_csv(RAP_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/RAP_final'))  

