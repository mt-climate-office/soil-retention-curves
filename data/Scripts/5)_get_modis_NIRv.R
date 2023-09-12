#Get Modis data 
#Date Modified: 3/17/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script calculates NIRv from MODIS
# Addiontal code indluces Landsat, EVI, ET, and GPP


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


#select years of analysis 

years = 2016:2022

#Download NIRv data 
for(i in 1:length(years)){
  dataset <- ee$ImageCollection('MODIS/006/MCD43A4')
  
  start_date = paste0(years[i],'-01-01')
  end_date = paste0(years[i],'-12-31')
  
  red_band <- dataset$
    select('Nadir_Reflectance_Band1')$
    filterDate(start_date, end_date)$
    toBands()
  
  red_qa_band <- dataset$
    select('BRDF_Albedo_Band_Mandatory_Quality_Band1')$
    filterDate(start_date, end_date)$
    toBands()
  
  nir_band <- dataset$
    select('Nadir_Reflectance_Band2')$
    filterDate(start_date, end_date)$
    toBands()
  
  nir_qa_band <- dataset$
    select('BRDF_Albedo_Band_Mandatory_Quality_Band2')$
    filterDate(start_date, end_date)$
    toBands()
  
  #extract by points
  red = ee_extract(
    #MODIS multiband image
    x = red_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 500,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  red_qa = ee_extract(
    #MODIS multiband image
    x = red_qa_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 500,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  #extract by points
  nir = ee_extract(
    #MODIS multiband image
    x = nir_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 500,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  nir_qa = ee_extract(
    #MODIS multiband image
    x = nir_qa_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 500,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  red_clean = red %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'red_sr', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'),
           red_sr = red_sr*0.0001)
  
  red_qa_clean = red_qa %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'red_sr_qa', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'))
  
  nir_clean = nir %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'nir_sr', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'),
           nir_sr = nir_sr*0.0001)
  
  nir_qa_clean = nir_qa %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'nir_sr_qa', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'))
  
  joined = left_join(red_clean, nir_clean, by = c('station',  'name', 'date_installed', 'sub_network', 'elevation', 'mesowest_id', 'time')) %>%
    left_join(., red_qa_clean, by = c('station',  'name', 'date_installed', 'sub_network', 'elevation', 'mesowest_id', 'time')) %>%
    left_join(., nir_qa_clean, by = c('station',  'name', 'date_installed', 'sub_network', 'elevation', 'mesowest_id', 'time')) %>%
    mutate(red_sr = ifelse(red_sr_qa == 0, red_sr, NA),
           nir_sr = ifelse(nir_sr_qa == 0, nir_sr, NA),
           ndvi = (nir_sr - red_sr)/(nir_sr + red_sr),
           ndvi = ifelse(ndvi > 1, 1, ndvi),
           NIRv = (ndvi - 0.08) * nir_sr,
           NIRv_rolling = zoo::rollmean(NIRv, 7, align = 'center', na.rm = T, fill = NA),
           ndvi_rolling = zoo::rollmean(ndvi, 7, align = 'center', na.rm = T, fill = NA))
  
  write_csv(joined, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/NIRv/NIRv_',years[i],'.csv'))
}


binded = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/NIRv', full.names = T) %>%
  purrr::map(., read_csv) %>%
  bind_rows()

test = joined %>%
  filter(station == 'blmroyno')
plot(test$NIRv)
lines(test$NIRv_rolling, add = T)




# Landsat data

for(i in 1:length(years)){
  dataset <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')
  
  start_date = paste0(years[i],'-04-01')
  end_date = paste0(years[i],'-04-30')
  
  red_band <- dataset$
    select('SR_B4')$
    filterDate(start_date, end_date)$
    toBands()
  
  nir_band <- dataset$
    select('SR_B5')$
    filterDate(start_date, end_date)$
    toBands()
  
  #extract by points
  red = ee_extract(
    #MODIS multiband image
    x = red_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 30,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  #extract by points
  nir = ee_extract(
    #MODIS multiband image
    x = nir_band,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (250m = native)
    scale = 30,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  red_clean = red %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'red_sr', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'),
           red_sr = red_sr*0.0001)
  
  nir_clean = nir %>%
    st_drop_geometry() %>%
    dplyr::select(-gwic_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name, date_installed, sub_network, elevation, mesowest_id), values_to = 'nir_sr', names_to = 'time') %>%
    mutate(time = as.Date(substr(time, 2, 11), format = '%Y_%m_%d'),
           nir_sr = nir_sr*0.0001)

  joined = left_join(red_clean, nir_clean, by = c('station',  'name', 'date_installed', 'sub_network', 'elevation', 'mesowest_id', 'time')) %>%
    mutate(red_sr = ifelse(red_sr_qa == 0, red_sr, NA),
           nir_sr = ifelse(nir_sr_qa == 0, nir_sr, NA),
           ndvi = (nir_sr - red_sr)/(nir_sr + red_sr),
           ndvi = ifelse(ndvi > 1, 1, ndvi),
           NIRv = (ndvi - 0.08) * nir_sr,
           NIRv_rolling = zoo::rollmean(NIRv, 7, align = 'center', na.rm = T, fill = NA),
           ndvi_rolling = zoo::rollmean(ndvi, 7, align = 'center', na.rm = T, fill = NA))
  
  write_csv(joined, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/NIRv/NIRv_Landsat_',years[i],'.csv'))
}


binded = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/NIRv', full.names = T) %>%
  purrr::map(., read_csv) %>%
  bind_rows()

test = joined %>%
  filter(station == 'blmroyno')
plot(test$NIRv)
lines(test$NIRv_rolling, add = T)





#EVI Data

years = 2017:2022

for(i in 1:length(years)){
  
  start_date = paste0(years[i],'-01-01')
  end_date = paste0(years[i],'-12-31')


evi = ee$ImageCollection("MODIS/MOD09GA_006_EVI")$
  filter(ee$Filter$date(start_date, end_date))$
  filter(ee$Filter$dayOfYear(152, 273))$
  toBands()

ET_band <- dataset$
  select('ET')$
  filterDate(start_date, end_date)$
  toBands()


#extract by points
point_extract = ee_extract(
  #ee image
  #  x = evi,
  
  #RAP image
  x = evi,
  #sf object (no need to convert to ee object now)
  y = mesonet,
  #sampling function - mean works, you are extracting at native res
  fun = ee$Reducer$mean(),
  #define resolution for extraction (30m = native)
  scale = 30,
  #tell rgee the shp file is a sf object
  sf = T
) 


evi_data = point_extract %>%
  st_drop_geometry() %>%
  select(-c(date_installed, sub_network, elevation, mesowest_id, gwic_id)) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(station, name), names_to = 'date') %>%
  mutate(date = substr(date, 2, 11) %>%
           as.Date(., format = '%Y_%m_%d'))


write_csv(evi_data, paste0('/home/kjamerson/pedotransfer/Research/EVI/evi_data_',years[i],'.csv'))

}


# MODIS ET
years = 2016:2022
years = 2016
i = 1
for(i in 1:length(years)){
  
  dataset <- ee$ImageCollection('MODIS/006/MOD16A2')
  
  start_date = paste0(years[i],'-01-01')
  end_date = paste0(years[i],'-12-31')
  
  
 ET <- dataset$
    select('ET')$
    filterDate(start_date, end_date)$
    toBands()
  
  #extract by points
  point_extract = ee_extract(
    #ee image
    #  x = evi,
    
    #RAP image
    x = ET,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (30m = native)
    scale = 30,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  ET_data = point_extract %>%
    st_drop_geometry() %>%
    select(-c(date_installed, sub_network, elevation, mesowest_id, gwic_id)) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name), names_to = 'date') %>%
    mutate(date = substr(date, 2, 11) %>%
             as.Date(., format = '%Y_%m_%d'))
  
  
  write_csv(ET_data, paste0('/home/kjamerson/pedotransfer/Research/ET/ET_data_',years[i],'.csv'))
  
}


# MODIS GPP
years = 2016:2022
years = 2016
i = 1
for(i in 1:length(years)){
  
  dataset <- ee$ImageCollection('MODIS/006/MOD17A2H')
  
  start_date = paste0(years[i],'-01-01')
  end_date = paste0(years[i],'-12-31')
  
  
  GPP <- dataset$
    select('Gpp')$
    filterDate(start_date, end_date)$
    toBands()
  
  #extract by points
  point_extract = ee_extract(
    #ee image
    #  x = evi,
    
    #RAP image
    x = GPP,
    #sf object (no need to convert to ee object now)
    y = mesonet,
    #sampling function - mean works, you are extracting at native res
    fun = ee$Reducer$mean(),
    #define resolution for extraction (30m = native)
    scale = 30,
    #tell rgee the shp file is a sf object
    sf = T
  ) 
  
  GPP_data = point_extract %>%
    st_drop_geometry() %>%
    select(-c(date_installed, sub_network, elevation, mesowest_id, gwic_id)) %>%
    as_tibble() %>%
    pivot_longer(cols = -c(station, name), names_to = 'date') %>%
    mutate(date = substr(date, 2, 11) %>%
             as.Date(., format = '%Y_%m_%d'))
  
  
  write_csv(GPP_data, paste0('/home/kjamerson/pedotransfer/Research/GPP/GPP_data_',years[i],'.csv'))
  
}
