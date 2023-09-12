#Convert VWC to psi 
#Date Modified: 3/16/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script converts vwc measurements from the Montana Mesonet weather stations into soil water potential measurements. 
#This scripts reads the RDS files computed in Script 1. 
#This script is broken up into 2 inch, 4 inch, 8 inch, 20 inch and and computes data one depth at a time and stores them in separate folders. 
#Script #3 will combine them into one folder. 
#There is also code to compute the 36 inch depths but script #3 only uses shallow (2-4 inches), mid (8 inches), and deep (20 inches). 


#Load libraries 
library(tidyverse)
library(data.table)
library(magrittr)


# find all data
files_short = list.files('/home/kjamerson/pedotransfer/Research/SWRC_Plots', full.names = F, pattern = ".RDS$")
files_full = list.files('/home/kjamerson/pedotransfer/Research/SWRC_Plots', full.names = T, pattern = ".RDS$")
  
# filter models for depths
depths = substr(files_short, 9, 10)

#compute valid stations
stations = substr(files_short, 1, 8) %>%
  unique()

#Some of the shallow depths are recording at 2 inches and some are at 4 inches. Index 2 inches first, then 4 inches
#then manually add them to the same folder


index_depths = which(depths == '02')


# filter for valid depth indicies
files_short_filtered = files_short[index_depths]

files_full_filtered = files_full[index_depths]

stations = substr(files_short_filtered, 1, 8) %>%
  unique()


#Write directory for 2 inch depth

write_csv(stations %>% as_tibble(), file = '/home/kjamerson/pedotransfer/Research/Water_Potential/2_Inch/valid_station.csv')


# filter valid mesonet data based on what stations have been processed 
valid_vwc_data = list.files('/home/kjamerson/pedotransfer/Research/mesonet', full.names = T) %>%
  as_tibble() %>%
  filter(value %like% paste0(stations, collapse = "|")) %$%
  value


for(i in 1:length(stations)){

 # import soil moisture and filter for 2 in depth
  temp_vwc = read_csv(valid_vwc_data[i]) %>%
     filter(name %in% c('soilwc04'))
  
  # bring in soil temperature data 
  temp_soilt_04 = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilt_04')) %>%
    select (station_key, datetime, value ) 
  colnames(temp_soilt_04) = c("station_key", "datetime", "soilt_04_C")
  
  
# Calculate VPD
  #temp units = C, rh units = percent
  temp_rh_t = (left_join(read_csv(valid_vwc_data[i]) %>%
                        filter(name %in% c('rel_humi')) %>%
                        select(station_key, datetime, name, value),
                      read_csv(valid_vwc_data[i]) %>%
                        filter(name %in% c('air_temp')) %>%
                        select(station_key, datetime, name, value),
                      by = c('station_key', 'datetime'))) %>%
    mutate(air_temp = value.y, rel_humi = value.x)%>%
    mutate(rel_humi = ifelse(rel_humi>1, rel_humi/100, rel_humi)) %>%
    mutate(sat_VP = 0.6108*exp((17.27*air_temp)/(air_temp + 237.3)),
           vpd_kPa = (1-rel_humi)*sat_VP) %>%     
    select(station_key, datetime, vpd_kPa)
 
  
  #import water retention model
  temp_model = readRDS(files_full_filtered[i])
  
  #define coeffitients 
  model_coef = coef(temp_model[[1]])


#Convert vwc to psi, based on best model 
  
  if(temp_model[[5]] == "Fredlund-Xing Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((model_coef['h'])*(((exp(1)^(((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1/model_coef['m'])) - (exp(1)))^(1/model_coef['n']))))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize columns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
    
    #left join VPD here with mutate
    left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data 
   final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    #write our data
   
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/2_Inch/water_potential_', stations[i], '04.csv'))
  }

  
  if(temp_model[[5]] == "Van Genechten Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = (((((model_coef['s']-model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1-(1/model_coef['n']))) - 1)/model_coef['a'])^(1/model_coef['n'])) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
   
    #Left join with soil temp data  
    final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    #write our data
    
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/2_Inch/water_potential_', stations[i], '04.csv'))
  }
  
  if(temp_model[[5]] == "Brooks Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = (model_coef['a']*((temp_vwc$value - model_coef['r'])/(model_coef['s']-model_coef['r']))^(-1/model_coef['n']))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data 
    final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    #write our data
    
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/2_Inch/water_potential_', stations[i], '04.csv'))
  }


  if(temp_model[[5]] == "Kosugi Model"){
    
    which_res <- function(data) {
      model_data <- temp_model[[3]]
      kPa = rep(NA, length (temp_vwc$value))
      
      for(i in 1:length(temp_vwc$value)){ 
        if(is.na(temp_vwc$value[i])) {}
        else {
          
          res = abs((temp_vwc$value[i] - model_data$fit_VWC))
          kPa[i] = model_data$kPa[which(res == min(res))]
          
        }
      }
      new_data = tibble(psi_kPa = kPa) %>%
     mutate(
       datetime = temp_vwc$datetime,
       station_key = temp_vwc$station_key,
       vwc = temp_vwc$value) %>%
        #re organize collumns
        select(station_key, datetime, psi_kPa, vwc) %>%
        #filter data for max obs kPa value
        filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
        
        #  left join VPD here with mutate
        
        left_join(temp_rh_t, new_data, by = c("station_key" = 'station_key', "datetime"="datetime")) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa)
      
      #left join with soil temp data
      soil_data = left_join (new_data,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
      
      return(soil_data)
    }
    
last_step <- which_res()
    #write our data
    write_csv(last_step, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/2_Inch/water_potential_', stations[i], '04.csv'))
    
    
  }
  print(temp_model[[4]])
}


#############################################################################
# filter out 4 inch depths
index_depths = which(depths == '04')

files_short_filtered = files_short[index_depths]

files_full_filtered = files_full[index_depths]

stations = substr(files_short_filtered, 1, 8) %>%
  unique()

write_csv(stations %>% as_tibble(), file = '/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch/valid_station.csv')



valid_vwc_data = list.files('/home/kjamerson/pedotransfer/Research/mesonet', full.names = T) %>%
  as_tibble() %>%
  filter(value %like% paste0(stations, collapse = "|")) %$%
  value


for(i in 1:length(stations)){
  # import soil moisture and filter for 4 in depth
  temp_vwc = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilwc04'))

  # bring in soil temperature data 
  
  temp_soilt_04 = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilt_04')) %>%
    select (station_key, datetime, value ) 
  colnames(temp_soilt_04) = c("station_key", "datetime", "soilt_04_C")
  
  #Compute VPD
  #temp units = C, rh units = percent
  temp_rh_t = (left_join(read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('rel_humi')) %>%
                           select(station_key, datetime, name, value),
                         read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('air_temp')) %>%
                           select(station_key, datetime, name, value),
                         by = c('station_key', 'datetime'))) %>%
    mutate(air_temp = value.y, rel_humi = value.x)%>%
    mutate(rel_humi = ifelse(rel_humi>1, rel_humi/100, rel_humi)) %>%
    mutate(sat_VP = 0.6108*exp((17.27*air_temp)/(air_temp + 237.3)),
           vpd_kPa = (1-rel_humi)*sat_VP) %>%     
    select(station_key, datetime, vpd_kPa)
  
  
  #import water retention model
  temp_model = readRDS(files_full_filtered[i])
  
  #define coeffitients 
  model_coef = coef(temp_model[[1]])
  
  
  
  if(temp_model[[5]] == "Fredlund-Xing Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((model_coef['h'])*(((exp(1)^(((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1/model_coef['m'])) - (exp(1)))^(1/model_coef['n']))))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
    
    #left join VPD here with mutate
    left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    
    #write our data
    
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch/water_potential_', stations[i], '04.csv'))
 
  }
  
  
  if(temp_model[[5]] == "Van Genechten Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((((((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^((model_coef['n']/(model_coef['n']-1))))-1)^(1/model_coef['n']))/model_coef['a'])) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #  left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch/water_potential_', stations[i], '04.csv'))
  }
  
  if(temp_model[[5]] == "Brooks Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = (model_coef['a']*((temp_vwc$value - model_coef['r'])/(model_coef['s']-model_coef['r']))^(-1/model_coef['n']))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
    
    #write our data
    
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch/water_potential_', stations[i], '04.csv'))
  }
  

  if(temp_model[[5]] == "Kosugi Model"){
    
    which_res <- function(data) {
      model_data <- temp_model[[3]]
      kPa = rep(NA, length (temp_vwc$value))
      
      for(i in 1:length(temp_vwc$value)){ 
        if(is.na(temp_vwc$value[i])) {}
        else {
          
          res = abs((temp_vwc$value[i] - model_data$fit_VWC))
          kPa[i] = model_data$kPa[which(res == min(res))]
          
        }
      }
      new_data = tibble(psi_kPa = kPa) %>%
        mutate(
          datetime = temp_vwc$datetime,
          station_key = temp_vwc$station_key,
          vwc = temp_vwc$value) %>%
        #re organize collumns
        select(station_key, datetime, psi_kPa, vwc) %>%
        #filter data for max obs kPa value
        filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
        
        #left join VPD here with mutate
        left_join(temp_rh_t, new_data, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
      
      #left join with soil temp data
      soil_data = left_join (new_data,temp_soilt_04, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_04_C)
      
      return(soil_data)
    }
    
   last_step <- which_res()
    
    
    
    #write our data
    write_csv(last_step, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch/water_potential_', stations[i], '04.csv'))
    
    
  }
  print(temp_model[[4]])
  
}

########################################################################################################################
# Filter out 8 inch depths 

index_depths = which(depths == '08')

# filter for valid depth indices
files_short_filtered = files_short[index_depths]

files_full_filtered = files_full[index_depths]

stations = substr(files_short_filtered, 1, 8) %>%
  unique()

write_csv(stations %>% as_tibble(), file = '/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch/valid_station.csv')


valid_vwc_data = list.files('/home/kjamerson/pedotransfer/Research/mesonet', full.names = T) %>%
  as_tibble() %>%
  filter(value %like% paste0(stations, collapse = "|")) %$%
  value


for(i in 1:length(stations)){
  
 
  # import soil moisture and filter for 8 in depth
  temp_vwc = read_csv(valid_vwc_data[i]) %>%
  filter(name %in% c('soilwc08'))
  
  temp_soilt_08 = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilt_08')) %>%
    select (station_key, datetime, value ) 
  colnames(temp_soilt_08) = c("station_key", "datetime", "soilt_08_C")

  #Calculate VPD
  #temp units = C, rh units = percent
  temp_rh_t = (left_join(read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('rel_humi')) %>%
                           select(station_key, datetime, name, value),
                         read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('air_temp')) %>%
                           select(station_key, datetime, name, value),
                         by = c('station_key', 'datetime'))) %>%
    mutate(air_temp = value.y, rel_humi = value.x)%>%
    mutate(rel_humi = ifelse(rel_humi>1, rel_humi/100, rel_humi)) %>%

    mutate(sat_VP = 0.6108*exp((17.27*air_temp)/(air_temp + 237.3)),
           vpd_kPa = (1-rel_humi)*sat_VP) %>%     
    select(station_key, datetime, vpd_kPa)
  
  
  #import water retention model
  temp_model = readRDS(files_full_filtered[i])
  
  #define coeffitients 
  model_coef = coef(temp_model[[1]])
 
  
  if(temp_model[[5]] == "Fredlund-Xing Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((model_coef['h'])*(((exp(1)^(((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1/model_coef['m'])) - (exp(1)))^(1/model_coef['n']))))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_08, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_08_C)
    
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch/water_potential_', stations[i], '08.csv'))

  }
  
  
  if(temp_model[[5]] == "Van Genechten Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((((((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^((model_coef['n']/(model_coef['n']-1))))-1)^(1/model_coef['n']))/model_coef['a'])) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #Left join with soil temp data 
    final_data = left_join (temp_kPa,temp_soilt_08, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_08_C)
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch/water_potential_', stations[i], '08.csv'))
  }
  
  if(temp_model[[5]] == "Brooks Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = (model_coef['a']*((temp_vwc$value - model_coef['r'])/(model_coef['s']-model_coef['r']))^(-1/model_coef['n']))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_08, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_08_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch/water_potential_', stations[i], '08.csv'))
  }
  
  
  if(temp_model[[5]] == "Kosugi Model"){
    
    which_res <- function(data) {
      model_data <- temp_model[[3]]
      kPa = rep(NA, length (temp_vwc$value))
      
      for(i in 1:length(temp_vwc$value)){ 
        if(is.na(temp_vwc$value[i])) {}
        else {
          
          res = abs((temp_vwc$value[i] - model_data$fit_VWC))
          kPa[i] = model_data$kPa[which(res == min(res))]
          
        }
      }
      new_data = tibble(psi_kPa = kPa) %>%
        mutate(
          datetime = temp_vwc$datetime,
          station_key = temp_vwc$station_key,
          vwc = temp_vwc$value) %>%
        #re organize collumns
        select(station_key, datetime, psi_kPa, vwc) %>%
        #filter data for max obs kPa value
        filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
        
        #left join VPD here with mutate
        left_join(temp_rh_t, new_data, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
      
      #left join with soil temp data
      soil_data = left_join (new_data,temp_soilt_08, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_08_C)
      
      return(soil_data)
    }
    
    last_step <- which_res()
    
    #write our data
    write_csv(last_step, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch/water_potential_', stations[i], '08.csv'))
  }
  print(temp_model[[4]])
}

###################################################################################################################################
#20 Inch Depth 
index_depths = which(depths == '20')

# filter for valid depth indicies
files_short_filtered = files_short[index_depths]

files_full_filtered = files_full[index_depths]

stations = substr(files_short_filtered, 1, 8) %>%
  unique()

write_csv(stations %>% as_tibble(), file = '/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch/valid_station.csv')

valid_vwc_data = list.files('/home/kjamerson/pedotransfer/Research/mesonet', full.names = T) %>%
  as_tibble() %>%
  filter(value %like% paste0(stations, collapse = "|")) %$%
  value


for(i in 1:length(stations)){
  
  # import soil moisture and filter for 20 in depth
  temp_vwc = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilwc20'))
  
  temp_soilt_20 = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilt_20')) %>%
    select (station_key, datetime, value ) 
  colnames(temp_soilt_20) = c("station_key", "datetime", "soilt_20_C")
 
  #Calculate VPD
  #temp units = C, rh units = percent
  temp_rh_t = (left_join(read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('rel_humi')) %>%
                           select(station_key, datetime, name, value),
                         read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('air_temp')) %>%
                           select(station_key, datetime, name, value),
                         by = c('station_key', 'datetime'))) %>%
    mutate(air_temp = value.y, rel_humi = value.x)%>%
    mutate(rel_humi = ifelse(rel_humi>1, rel_humi/100, rel_humi)) %>%
    mutate(sat_VP = 0.6108*exp((17.27*air_temp)/(air_temp + 237.3)),
           vpd_kPa = (1-rel_humi)*sat_VP) %>%     
    select(station_key, datetime, vpd_kPa)
  
  
  #import water retention model
  temp_model = readRDS(files_full_filtered[i])
  
  #define coeffitients 
  model_coef = coef(temp_model[[1]])
  

  if(temp_model[[5]] == "Fredlund-Xing Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((model_coef['h'])*(((exp(1)^(((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1/model_coef['m'])) - (exp(1)))^(1/model_coef['n']))))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_20, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_20_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch/water_potential_', stations[i], '20.csv'))

  }
  
  
  if(temp_model[[5]] == "Van Genechten Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((((((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^((model_coef['n']/(model_coef['n']-1))))-1)^(1/model_coef['n']))/model_coef['a'])) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_20, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_20_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch/water_potential_', stations[i], '20.csv'))

  }
  
  if(temp_model[[5]] == "Brooks Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = (model_coef['a']*((temp_vwc$value - model_coef['r'])/(model_coef['s']-model_coef['r']))^(-1/model_coef['n']))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_20, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_20_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch/water_potential_', stations[i], '20.csv'))
    print(i)
  }
  
  if(temp_model[[5]] == "Kosugi Model"){
    
    which_res <- function(data) {
      model_data <- temp_model[[3]]
      kPa = rep(NA, length (temp_vwc$value))
      
      for(i in 1:length(temp_vwc$value)){ 
        if(is.na(temp_vwc$value[i])) {}
        else {
          
          res = abs((temp_vwc$value[i] - model_data$fit_VWC))
          kPa[i] = model_data$kPa[which(res == min(res))]
          
        }
      }
      new_data = tibble(psi_kPa = kPa) %>%
        mutate(
          datetime = temp_vwc$datetime,
          station_key = temp_vwc$station_key,
          vwc = temp_vwc$value) %>%
        #re organize collumns
        select(station_key, datetime, psi_kPa, vwc) %>%
        #filter data for max obs kPa value
        filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
        
        #left join VPD here with mutate
        left_join(temp_rh_t, new_data, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
      
     #lelft join with soil temp data
       soil_data = left_join (new_data,temp_soilt_20, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_20_C)
      
      return(soil_data)
    }
    
    last_step <- which_res()
  
    #write our data
    write_csv(last_step, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch/water_potential_', stations[i], '20.csv'))
  }
 
  }
  print(temp_model[[4]])
  
  


##################################################################################################################################
  # 36 inch depth
index_depths = which(depths == '36')

# filter for valid depth indicies
files_short_filtered = files_short[index_depths]

files_full_filtered = files_full[index_depths]

stations = substr(files_short_filtered, 1, 8) %>%
  unique()

write_csv(stations %>% as_tibble(), file = '/home/kjamerson/pedotransfer/Research/Water_Potential/36_Inch/valid_station.csv')

valid_vwc_data = list.files('/home/kjamerson/pedotransfer/Research/mesonet', full.names = T) %>%
  as_tibble() %>%
  filter(value %like% paste0(stations, collapse = "|")) %$%
  value

for(i in 1:length(stations)){
  # import soil moisture and filter for 36 in depth
    temp_vwc = read_csv(valid_vwc_data[i]) %>%
    filter(name %in% c('soilwc36'))
    
    temp_soilt_36 = read_csv(valid_vwc_data[i]) %>%
      filter(name %in% c('soilt_36')) %>%
      select (station_key, datetime, value ) 
    colnames(temp_soilt_36) = c("station_key", "datetime", "soilt_36_C")
  
  #Calculate VPD
  #temp units = C, rh units = percent
  temp_rh_t = (left_join(read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('rel_humi')) %>%
                           select(station_key, datetime, name, value),
                         read_csv(valid_vwc_data[i]) %>%
                           filter(name %in% c('air_temp')) %>%
                           select(station_key, datetime, name, value),
                         by = c('station_key', 'datetime'))) %>%
    mutate(air_temp = value.y, rel_humi = value.x)%>%
    mutate(rel_humi = ifelse(rel_humi>1, rel_humi/100, rel_humi)) %>%
    #If using bigleaf function, do not divide by 100
    #mutate(vpd = rH.to.VPD(rel_humi, air_temp)) %>% #Bigleaf function that calculates VPD
    #mutate (vpd = VPD*(-1)) %>% #Negative values, need to be postive? 
    #Bigleaf function that calculates VPD
    mutate(sat_VP = 0.6108*exp((17.27*air_temp)/(air_temp + 237.3)),
           vpd_kPa = (1-rel_humi)*sat_VP) %>%     
    select(station_key, datetime, vpd_kPa)
  
  #import water retention model
  temp_model = readRDS(files_full_filtered[i])
  
  #define coeffitients 
  model_coef = coef(temp_model[[1]])
  
  if(temp_model[[5]] == "Fredlund-Xing Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((model_coef['h'])*((exp(1)^(((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^(1/model_coef['m'])) - (exp(1))^(1/model_coef['n']))))) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize collumns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_36, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_36_C)
    
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/36_Inch/water_potential_', stations[i], '36.csv'))
  }
  
  
  if(temp_model[[5]] == "Van Genechten Model"){
    #compute inverse model
    temp_kPa = tibble(psi_kPa = ((((((model_coef['s'] - model_coef['r'])/(temp_vwc$value - model_coef['r']))^((model_coef['n']/(model_coef['n']-1))))-1)^(1/model_coef['n']))/model_coef['a'])) %>%
      #add in data from mesonet file
      mutate(datetime = temp_vwc$datetime,
             station_key = temp_vwc$station_key,
             vwc = temp_vwc$value) %>%
      #re organize columns
      select(station_key, datetime, psi_kPa, vwc) %>%
      #filter data for max obs kPa value
      filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
      
      #left join VPD here with mutate
      left_join(temp_rh_t, temp_kPa, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
    
    #left join with soil temp data
    final_data = left_join (temp_kPa,temp_soilt_36, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
      select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_36_C)
    
    #write our data
    write_csv(final_data, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/36_Inch/water_potential_', stations[i], '36.csv'))

  }
  
  if(temp_model[[5]] == "Kosugi Model"){
    
    which_res <- function(data) {
      model_data <- temp_model[[3]]
      kPa = rep(NA, length (temp_vwc$value))
      
      for(i in 1:length(temp_vwc$value)){ 
        if(is.na(temp_vwc$value[i])) {}
        else {
          
          res = abs((temp_vwc$value[i] - model_data$fit_VWC))
          kPa[i] = model_data$kPa[which(res == min(res))]
          
        }
      }
      new_data = tibble(psi_kPa = kPa) %>%
        mutate(
          datetime = temp_vwc$datetime,
          station_key = temp_vwc$station_key,
          vwc = temp_vwc$value) %>%
        #re organize collumns
        select(station_key, datetime, psi_kPa, vwc) %>%
        #filter data for max obs kPa value
        filter(psi_kPa < max(temp_model[[2]]$raw_kPa)) %>%
        
        #left join VPD here with mutate
        left_join(temp_rh_t, new_data, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa) 
      
      #left join with soil temp data
      soil_data = left_join (new_data,temp_soilt_36, by = c('station_key' = 'station_key', 'datetime' = 'datetime')) %>%
        select(station_key, datetime, psi_kPa, vwc, vpd_kPa, soilt_36_C)
      
      return(soil_data)
    }
    
    last_step <- which_res()
    
    
    #write our data
    write_csv(last_step, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/36_Inch/water_potential_', stations[i], '36.csv'))
  }
  
  }
  print(i) 
  

