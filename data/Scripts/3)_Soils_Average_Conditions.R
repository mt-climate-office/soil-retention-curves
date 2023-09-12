#Average soil conditions 
#Date Modified: 3/17/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script calculates daily averages for each station of soil psi, VPD, vwc, and soil temp for shallow (2-4 inches), mid (8 inches) and deep (20 inches) soil depths.  


library(matrixStats)
library(readr)
library(tidyverse)



# find all data
Depth04 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch', full.names = T)
Depth08 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch', full.names = T)
Depth20 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch', full.names = T)

#Find sites that don't have all 3 depths and remove
names_4 = substr(Depth04, 78, 85)
names_8 = substr(Depth08, 78, 85)
names_20 = substr(Depth20, 79, 86)


names4_names8 = which(! names_4 %in% names_8)

delete_from_4 = names_4[names4_names8]

names8_names20 = which(! names_8 %in% names_20)

delete_from_8 = names_8[names8_names20]


#join all depths together in one file
files = tibble(files_04 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/4_Inch', full.names = T),
              files_08 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/8_Inch', full.names = T),
              files_20 = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/20_Inch', full.names = T)) %>%
  mutate(stations_4 = substr(files_04, 78, 85),
         stations_8 = substr(files_08, 78, 85),
         stations_20 = substr(files_20, 79, 86))


#make sure that all depths have the same number of stations 
if(all.equal(files$stations_4,files$stations_8) == T){
  if(all.equal(files$stations_8,files$stations_20) == T){
    if(all.equal(files$stations_4,files$stations_20) == T){
      print('All equal')
    }
  }
}


#Create station list
  stations = files$stations_4




for(i in 1:length(stations)) {
 
  #Compute daily averages for shallow depths (4 Inch) 
  Daily_Data_Sh <- read_csv(files$files_04[i]) %>%
    mutate(date = as.Date(datetime)) %>%
    pivot_longer(cols = -c(station_key, datetime, date)) %>%
    #group_by date and name (variable)
    group_by(date, name) %>%
    #summarize to daily
    summarise(daily_value = mean(value)) %>%
    #pivot back to wide form
    pivot_wider(names_from = name, values_from = daily_value) %>%
    #add back in station ID
    mutate(station_key = stations[i]) %>%
    #rearrange columns
    select(station_key, date, psi_kPa, vpd_kPa, vwc, soilt_04_C)
  colnames(Daily_Data_Sh) = c("station", "date", "psi_kPa_sh","vpd_kPa_sh",'vwc_sh',"soilt_C_sh") 
  
 #compute daily averages for mid depths (8 Inch) 
  Daily_Data_Mid <- read_csv(files$files_08[i]) %>%
    mutate(date = as.Date(datetime)) %>%
    pivot_longer(cols = -c(station_key, datetime, date)) %>%
    #group_by date and name (variable)
    group_by(date, name) %>%
    #summarize to daily
    summarise(daily_value = mean(value)) %>%
    #pivot back to wide form
    pivot_wider(names_from = name, values_from = daily_value) %>%
    #add back in station ID
    mutate(station_key = stations[i]) %>%
    #rearange collumns
    select(station_key, date, psi_kPa, vpd_kPa, vwc, soilt_08_C)
  colnames(Daily_Data_Mid) = c("station", "date", "psi_kPa_mid","vpd_kPa_mid",'vwc_mid',"soilt_C_mid")

  #Join shallow and mid depths 
  Shallow_Mid = left_join(Daily_Data_Sh, Daily_Data_Mid, by = c('station' = 'station', 'date' = 'date')) 
  
#Compute daily averages for deep depths (20 inch)
  Daily_Data_Deep <- read_csv(files$files_20[i]) %>%
    mutate(date = as.Date(datetime)) %>%
    pivot_longer(cols = -c(station_key, datetime, date)) %>%
    #group_by date and name (variable)
    group_by(date, name) %>%
    #summarize to daily
    summarise(daily_value = mean(value)) %>%
    #pivot back to wide form
    pivot_wider(names_from = name, values_from = daily_value) %>%
    #add back in station ID
    mutate(station_key = stations[i]) %>%
    #rearrange columns
    select(station_key, date, psi_kPa, vpd_kPa, vwc, soilt_20_C)
  colnames(Daily_Data_Deep) = c("station", "date", "psi_kPa_deep","vpd_kPa_deep",'vwc_deep',"soilt_C_deep")
  
  #Left join all depths and then compute means/ranges 
  All_Depths= left_join (Shallow_Mid, Daily_Data_Deep, by = c('station' = 'station', 'date' = 'date')) %>%
    mutate(psi_kPa_mean= ((psi_kPa_sh + psi_kPa_mid + psi_kPa_deep)/3)) %>%
    mutate(vpd_kPa_mean =((vpd_kPa_sh + vpd_kPa_mid + vpd_kPa_deep)/3)) %>%
    mutate(vwc_mean =((vwc_sh + vwc_mid + vwc_deep)/3)) %>%
    mutate(soilt_C_mean = ((soilt_C_sh + soilt_C_mid + soilt_C_deep)/3)) %>%
    mutate(psi_range = (pmax(psi_kPa_sh,psi_kPa_mid, psi_kPa_deep)) - pmin(psi_kPa_sh,psi_kPa_mid, psi_kPa_deep) ) %>%
    mutate(vpd_range = (pmax(vpd_kPa_sh,vpd_kPa_mid, vpd_kPa_deep)) - pmin(vpd_kPa_sh,vpd_kPa_mid, vpd_kPa_deep) ) %>%
    mutate(vwc_range = (pmax(vwc_sh, vwc_mid, vwc_deep)) - pmin(vwc_sh,vwc_mid, vwc_deep) ) %>%
    mutate(soilt_C_range = (pmax(soilt_C_sh, soilt_C_mid, soilt_C_deep)) - pmin(soilt_C_sh,soilt_C_mid, soilt_C_deep) ) 
  
  Final = na.omit(All_Depths)
  
  #Save as CSV file
  write_csv(Final, paste0('/home/kjamerson/pedotransfer/Research/Water_Potential/Avg_Soil/avg_conditions_',stations[i]))
  
  print(i)
} 


