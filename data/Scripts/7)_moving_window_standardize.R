#Correlation plots of NIRv by vegetation cover
#Date Modified: 3/20/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script calculated anomalies using a 31-day moving window analysis. The response variable is NIRv and the predictor variables are soil water potential, VWC, VPD, and soil temperature. 
#Calculates anomalies, runs a full GAM model, calculated PE anomaly plots, partial deviance change tables, and anomaly plots

#Download libraries 
library(tidyr)
library(patchwork)
library(magrittr)
library(ggplot2)
library(data.table)
library(sf)
library(dplyr)
library(ggspatial)
library(ggrepel)
library(DescTools)
library(readr)
library(knitr)
library(mgcv)
library(gratia)
library(tidyverse)

source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

#minimum data threshold
min_data_thresh = 31*5

moving_window_standardize = function(date, variable, min_data_thresh, dist = 'gam'){ # dist = 'glo' glo can handle negatives
  #define tibble for processing
  x = tibble(date = date,
             variable = variable) %>%
    #compute some date meta
    dplyr::mutate(yday = lubridate::yday(date),
                  year = lubridate::year(date)) %>%
    #arange by time
    arrange(date)
  
  #compute unique ydays
  ydays = unique(x$yday)
  
  #compute last of each unique yday as indicies
  indicies = vector()
  for(i in 1:length(ydays)){
    indicies[i] = which(x$yday == ydays[i]) %>% 
      max()
  }
  
  out = list()
  for(index in indicies){
    date_of_interest = x$date[index]
    #compute index specific window
    top_window = x$yday[index] + 15
    bottom_window = x$yday[index] - 15
    #correct for yday breaks ar 1 and 365
    if(top_window > 365){
      top_window = top_window - 365
      top_range = c(x$yday[index]:365,1:top_window)
    } else {
      top_range = x$yday[index]:top_window
    }
    if(bottom_window < 1){
      bottom_window = bottom_window + 365
      bottom_range = c(bottom_window:365, 1:x$yday[index])
    } else {
      bottom_range = bottom_window:x$yday[index]
    }
    #compute range of y days
    range = c(bottom_range, top_range) %>% unique()
    
    #define function to return ecdf of vector
    ecdf_anomaly = function(x){
      ecdf_x = ecdf(x) 
      return(ecdf_x(x))
    }
    
    #filter data for index
    if(dist == 'gam'){
      standard = x %>%
        filter(yday %in% range) %>%
        #return all data (static reference frame) and long cliamtology length
        mutate(drought_anomaly = gamma_fit_spi(variable, return_latest = F, climatology_length = Inf, export_opts = 'CDF'),
               ecdf_anomaly = ecdf_anomaly(variable))
    }
    
    if(dist == 'glo'){
      standard = x %>% 
        filter(yday %in% range) %>%
        #return all data (static reference frame) and long cliamtology length
        mutate(drought_anomaly = glo_fit_spei(variable, return_latest = F, climatology_length = Inf, export_opts = 'CDF'),
               ecdf_anomaly = ecdf_anomaly(variable))
    }
    
    #find index of centroid date (date of interest)
    if(length(standard$date) >= min_data_thresh){
      out[[which(index == indicies)]] = standard[which(standard$yday == lubridate::yday(date_of_interest)),]
    } else {
      out[[which(index == indicies)]] = standard[which(standard$yday == lubridate::yday(date_of_interest)),] %>%
        mutate(drought_anomaly = NA,
               ecdf_anomaly = NA)
    }
  }
  #if its the last index
  if(index == indicies[length(indicies)]){
    #function to remove empty elements from list
    is_na <- function(x) {
      anyNA(x)
    }
    #bind the final data frame together
    out_final = out %>% 
      #purrr::discard(., is_na) %>%
      bind_rows() %>%
      dplyr::arrange(date)
  }
  return(out_final)
}

  #Create nontin to filter out years where environmental stations were first installed 
`%notin%` = Negate(`%in%`)

#Read in mesonet data
mesonet = read_csv('https://mesonet.climate.umt.edu/api/v2/stations/?type=csv&clean=true') %>%
  mutate(start_year = lubridate::year(date_installed),
         name_start_year = paste0(station,'_', start_year))

#Read in MSE data
MSE_data = read_csv('/home/kjamerson/pedotransfer/Research/MSE_Model_Data')

#import NIRv data
NIRv_binded = list.files('/home/kjamerson/pedotransfer/Research/Water_Potential/NIRv', full.names = T) %>%
  purrr::map(., read_csv) %>%
  bind_rows() %>%
  drop_na() %>%
  select(station, time, NIRv) %>%
  mutate( year = lubridate::year(time)) %>%
  filter(year %in% c(2017:2021),
         NIRv > 0)
colnames(NIRv_binded) = c("station", "date", 'NIRv','year')


#Import soil data
soils_binded = list.files(path = '~/pedotransfer/Research/Water_Potential/Avg_Soil', full.names = T) %>%
  purrr::map(.,read_csv) %>%
  bind_rows() %>%
  drop_na() %>%
  select(station, date, psi_kPa_mid, vpd_kPa_mid, vwc_mid, soilt_C_mid)%>%
  mutate(year = lubridate::year(date))

NIRv_Data = soils_binded %>%
  left_join(  NIRv_binded, by = c('station' = 'station', 'date' = 'date', 'year'='year')) %>%
  drop_na()


#Import RAP Data
RAP = read.csv('/home/kjamerson/pedotransfer/Research/Water_Potential/RAP_final')

#Create groups of dominant vegetation types
RAP_Data = left_join((RAP %>%
                        arrange(year, desc(value)) %>%
                        group_by(station, year) %>%
                        slice(1)),
                     (RAP %>%
                        arrange(year, desc(value)) %>%
                        group_by(station, year) %>%
                        slice(2)), by = c('station', 'name', 'year')) %>%
  mutate(top_groups = paste0(land_cover.x, ', ' ,land_cover.y),
         sum_area = sum(value.x, value.y)) %>%
  select(station, name, year, top_groups, sum_area) %>%
  ungroup() %>%
  mutate(top_groups = ifelse(top_groups == 'PFG, AFG', 'PFG, AFG',
                             ifelse(top_groups == 'PFG, LTR', 'PFG, LTR',
                                    ifelse(top_groups == 'PFG, BGR', 'PFG, BGR',
                                           ifelse(top_groups == 'PFG, SHR', 'PFG, SHR', 'Other')))),
         top_groups = ifelse(sum_area < 50, 'No Dominant Class', top_groups))



lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}



valid_sites = MSE_data %>%
  filter(MSE_data$Depth == '08',
         MSE < 0.00015)



#Summarize weekly data for all vegetation types, filtered for the growing season
NIRv_Final = NIRv_Data %>%
  mutate(year = lubridate::year(date)) %>%
  filter(station %in% valid_sites$station) %>%
  group_by(station) %>%
  mutate(min_year = min(year)) %>%
  filter(min_year <= 2018) %>%
  ungroup()


min_data_thresh_relaxed = 3*31

#Calculate moving window analysis 
NIRv_Data_standardized = NIRv_Final %>%
  group_by(station) %>%
  mutate(vwc_mid_standard = moving_window_standardize(date = date,
                                                      variable = vwc_mid,
                                                      min_data_thresh = min_data_thresh_relaxed,
                                                      dist = 'gam')$drought_anomaly,
         vpd_kPa_mid_standard = moving_window_standardize(date = date,
                                                      variable = vpd_kPa_mid,
                                                      min_data_thresh = min_data_thresh_relaxed,
                                                      dist = 'gam')$drought_anomaly,
         psi_kPa_mid_standard = moving_window_standardize(date = date,
                                                      variable = psi_kPa_mid,
                                                      min_data_thresh = min_data_thresh_relaxed,
                                                      dist = 'gam')$drought_anomaly,
         soilt_C_mid_standard = moving_window_standardize(date = date,
                                                          variable = soilt_C_mid,
                                                          min_data_thresh = min_data_thresh_relaxed,
                                                          dist = 'glo')$drought_anomaly,
         NIRv_standard = moving_window_standardize(date = date,
                                                          variable = NIRv,
                                                          min_data_thresh = min_data_thresh_relaxed,
                                                          dist = 'gam')$drought_anomaly)



NIRv_Data_standardized_filtered = NIRv_Data_standardized %>%
  mutate(vwc_mid_standard = qnorm(vwc_mid_standard),
         vpd_kPa_mid_standard = qnorm(vpd_kPa_mid_standard),
         psi_kPa_mid_standard = qnorm(psi_kPa_mid_standard),
         soilt_C_mid_standard = qnorm(soilt_C_mid_standard),
         NIRv_standard = qnorm(NIRv_standard)) %>%
  mutate(month = lubridate::month(date))%>%
  filter(month %in% 5:9) 

k_ = 6


NIRv_full <- gam(NIRv_standard ~ s(psi_kPa_mid_standard, k = k_) + s(vpd_kPa_mid_standard, k = k_) + s(soilt_C_mid_standard, k = k_), 
                 data = NIRv_Data_standardized_filtered, dispersion = 0.5)



#Create Partial Effect plots for each predictor variable 
NIRv_pe_psi = smooth_estimates(NIRv_full, smooth = "s(psi_kPa_mid_standard)")
NIRv_pe_vpd = smooth_estimates(NIRv_full, smooth = "s(vpd_kPa_mid_standard)")
NIRv_pe_soilt = smooth_estimates(NIRv_full, smooth = "s(soilt_C_mid_standard)")

water_psi = NIRv_Data_standardized_filtered %>%
ggplot(aes(y = vwc_mid_standard, x= psi_kPa_mid_standard))+
  geom_point(size = 2, shape = 21) + 
  geom_smooth(method = 'lm')+
  geom_abline(intercept = 0, slope = -1, color = 'red', linetype = 'dashed')+
  theme_bw(base_size = 20)+
  labs(y = "VWC (m³/m³) Anomaly",
       title = expression("Subsurface "),
       x = expression("Ψ (- kPa) Anomaly"), element_text( size = 30))

ggsave(water_psi, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/water_psi_anomaly.png'), width = 12, height = 8, units = 'in')


NIRv_pe_psi_plot = NIRv_pe_psi %>%
  add_confint() %>%
  ggplot(aes(y = est, x = psi_kPa_mid_standard))+
  theme_bw(base_size = 20)+
  geom_hline(yintercept=c(0), linetype = "dashed", color = "red", size = 1)+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "black") +
  geom_line(colour = "blue", size = 1) +
  labs(y = "Partial Effect",
       title = expression("Parametric Anomaly NIRv"),
       x = expression("Ψ (- kPa) Anomaly"), element_text( size = 20))

NIRv_pe_vpd_plot = NIRv_pe_vpd %>%
  add_confint() %>%
  ggplot(aes(y = est, x = vpd_kPa_mid_standard))+
  theme_bw(base_size = 20)+
  geom_hline(yintercept=c(0), linetype = "dashed", color = "red", size = 1)+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "black") +
  geom_line(color = "blue", size = 1) +
  labs(y = "Partial effect",
       x = expression("VPD (kPa) Anomaly"), element_text(size = 30))

NIRv_pe_soilt_plot = NIRv_pe_soilt %>%
  add_confint() %>%
  ggplot(aes(y = est, x = soilt_C_mid_standard))+
  theme_bw(base_size = 20)+
  geom_hline(yintercept=c(0), linetype = "dashed", color = "red", size = 1)+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "black") +
  geom_line(color = "blue", size = 1) +
  labs(y = "Partial effect",
       x = expression("Soil Temperature (C) Anomaly"), element_text(size = 30))


NIRv_PE_Plots = NIRv_pe_psi_plot + NIRv_pe_vpd_plot +  NIRv_pe_soilt_plot +  plot_layout(nrow = 3)
NIRv_PE_Plots
ggsave(NIRv_PE_Plots, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/NIRv_anomaly_PE_Plots_8_in.png'), width = 9, height = 14, units = 'in')


NIRv_null_psi_kPa_mid = gam(NIRv_standard ~  s(vpd_kPa_mid_standard)  + s(soilt_C_mid_standard), sp = NIRv_full$sp[-1], data = NIRv_Data_standardized_filtered, dispersion = 0.5)
NIRv_null_vpd_kPa_mid = gam(NIRv_standard ~ s(psi_kPa_mid_standard)  + s(soilt_C_mid_standard), sp = NIRv_full$sp[-2], data = NIRv_Data_standardized_filtered, dispersion = 0.5)
NIRv_null_soilt_C_mid = gam(NIRv_standard ~ s(psi_kPa_mid_standard) + s(vpd_kPa_mid_standard) , sp = NIRv_full$sp[-3], data = NIRv_Data_standardized_filtered, dispersion = 0.5)


# Subtract the null models from the full models. 
a_NIRv = round((deviance(NIRv_null_psi_kPa_mid)-deviance(NIRv_full)), digits = 5)
b_NIRv = round((deviance(NIRv_null_vpd_kPa_mid)-deviance(NIRv_full)), digits =5)
c_NIRv = round((deviance(NIRv_null_soilt_C_mid)-deviance(NIRv_full)), digits = 5)

#Create dataframe of partial Deviance Change
predictor_names <- c( "Ψ (kPa) Anomaly","VPD (kPa) Anomaly", "Soil Temperature (C) Anomaly")
pdc_NIRv = data.frame( variables = predictor_names, pdc = c(a_NIRv,b_NIRv,c_NIRv))
colnames(pdc_NIRv) = c("Variables", "Partial Deviance Change")
pdc_NIRv


NIRv_Data_standardized_filtered %>%
  select(psi_kPa_mid_standard, NIRv_standard) %>%
  drop_na()

psi_lm = lm(NIRv_standard ~ psi_kPa_mid_standard,
            data = NIRv_Data_standardized_filtered) %>% summary()

psi_p = psi_lm$coefficients[2,4]
if(psi_p < 0.05){
  psi_p = '< 0.05'
}


vpd_lm = lm(NIRv_standard ~ vpd_kPa_mid_standard,
            data = NIRv_Data_standardized_filtered) %>% summary()

vpd_p = vpd_lm$coefficients[2,4]

if(vpd_p < 0.05){
  vpd_p = '< 0.05'
}

soilt_lm = lm(NIRv_standard ~ soilt_C_mid_standard,
              data = NIRv_Data_standardized_filtered) %>% summary()

soilt_p = soilt_lm$coefficients[2,4]

if(soilt_p < 0.05){
  soilt_p = '< 0.05'
}


NIRv_psi_kPa = 
  ggplot(NIRv_Data_standardized_filtered, aes(x= psi_kPa_mid_standard, y = NIRv_standard))+
  labs(x = 'Ψ (- kPa) Anomaly', y = 'NIRv Anomaly')+
  geom_point(size = 2, shape = 21)+
  scale_color_discrete(name = NULL)+
  geom_text(aes(x = 2.5, y = 2.5, label = paste0('r = ',
                                                 NIRv_Data_standardized_filtered %>%
                                                   select(psi_kPa_mid_standard, NIRv_standard) %>%
                                                   drop_na() %$%
                                                   cor(psi_kPa_mid_standard, NIRv_standard) %>% round(., 3))), size = 5)+
  geom_text(aes(x = 2.5, y = 2.9, label = paste0('p = ',
                                                 psi_p)), size = 5)+
  theme_bw(base_size = 14)+
  geom_smooth(method = 'gam', formula = y ~ s(x, k = k_))+
  geom_smooth(method = 'lm', color = 'red')+
  #geom_smooth(method = 'lm')+
  ggtitle("Anomaly Analysis")+
  xlim(-3,3)+
  ylim(-3,3)+
  theme(plot.title = element_text( size = 27),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))+
  geom_hline(yintercept = 0, color = 'green')

NIRv_vpd_kPa = 
  ggplot(NIRv_Data_standardized_filtered, aes(x= vpd_kPa_mid_standard, y = NIRv_standard))+
  labs(x = 'VPD (kPa) Anomaly', y = 'NIRv Anomaly')+
  geom_point(size = 2, shape = 21)+
  scale_color_discrete(name = NULL)+
  theme_bw(base_size = 14)+
  geom_text(aes(x = 2.5, y = 2.5, label = paste0('r = ',
                                                 NIRv_Data_standardized_filtered %>%
                                                   select(vpd_kPa_mid_standard, NIRv_standard) %>%
                                                   drop_na() %$%
                                                   cor(vpd_kPa_mid_standard, NIRv_standard) %>% round(., 3))), size = 5)+
  geom_text(aes(x = 2.5, y = 2.9, label = paste0('p = ',
                                                 vpd_p)), size = 5)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_smooth(method = 'gam', formula = y ~ s(x))+
  geom_smooth(method = 'lm', color = 'red')+
  xlim(-3,3)+
  ylim(-3,3)+
  theme(plot.title = element_text( size = 27),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))+
  geom_hline(yintercept = 0, color = 'green')

NIRv_soilt_C = 
  ggplot(NIRv_Data_standardized_filtered, aes(x= soilt_C_mid_standard, y = NIRv_standard))+
  labs(x = 'Soil Temp (C) Anomaly', y = 'NIRv Anomaly')+
  geom_point(size = 2, shape = 21)+
  scale_color_discrete(name = NULL)+
  theme_bw(base_size = 14)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_text(aes(x = 2.5, y = 2.5, label = paste0('r = ',
                                                 NIRv_Data_standardized_filtered %>%
                                                   select(soilt_C_mid_standard, NIRv_standard) %>%
                                                   drop_na() %$%
                                                   cor(soilt_C_mid_standard, NIRv_standard) %>% round(., 3))), size = 5)+
  geom_text(aes(x = 2.5, y = 2.9, label = paste0('p = ',
                                                 soilt_p)), size = 5)+
  geom_smooth(method = 'gam', formula = y ~ s(x))+
  geom_smooth(method = 'lm', color = 'red')+
  xlim(-3,3)+
  ylim(-3,3)+
  theme(plot.title = element_text( size = 27),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))+
  geom_hline(yintercept = 0, color = 'green')


NIRv_Anomaly = NIRv_psi_kPa + NIRv_vpd_kPa + NIRv_soilt_C + plot_layout(nrow = 3)
NIRv_Anomaly
ggsave(NIRv_Anomaly, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/NIRv_anomaly_plots.png'), width = 9, height = 14, units = 'in')



#Create dataframe of gam full summary info

NIRv_dev = round(summary(NIRv_full)$dev.expl, digits = 4)*100


#Create dataframe of gam full summary info
NIRv_df = as.data.frame(summary(NIRv_full)$s.table)
attributes(NIRv_df)$row.names = "NULL"

#Join dataframes together and create a table 
NIRv_table = cbind(pdc_NIRv, NIRv_df) %>%
  arrange(desc(`Partial Deviance Change`))%>%
  kable(caption = "NIRv standardized anomaly GAM model, Deviance Explained = 26.74%")



NIRv_table %>%
  kableExtra::kable_styling() %>%
  kableExtra::row_spec(1, 1, extra_css = "border-bottom: 1px solid")%>%
  kableExtra::save_kable(file = "/home/kjamerson/pedotransfer/Research/Final_Figures/NIRv_anomaly_TABLE_mid.png")

Sys.setenv("OPENSSL_CONF"="/dev/null")

EVI_null_psi_kPa_mid = gam(EVI_standard ~  s(vpd_kPa_mid_standard)  + s(soilt_C_mid_standard), sp = EVI_full$sp[-1], data = EVI_Data_standardized_filtered, dispersion = 0.5)
EVI_null_vpd_kPa_mid = gam(EVI_standard ~ s(psi_kPa_mid_standard)  + s(soilt_C_mid_standard), sp = EVI_full$sp[-2], data = EVI_Data_standardized_filtered, dispersion = 0.5)
EVI_null_soilt_C_mid = gam(EVI_standard ~ s(psi_kPa_mid_standard) + s(vpd_kPa_mid_standard) , sp = EVI_full$sp[-3], data = EVI_Data_standardized_filtered, dispersion = 0.5)

a_EVI = round((deviance(EVI_null_psi_kPa_mid)-deviance(EVI_full)), digits = 5)
b_EVI = round((deviance(EVI_null_vpd_kPa_mid)-deviance(EVI_full)), digits =5)
c_EVI = round((deviance(EVI_null_soilt_C_mid)-deviance(EVI_full)), digits = 5)

EVI_pdc = data.frame( variables = predictor_names, EVI_pdc = c(a_EVI,b_EVI,c_EVI))
colnames(EVI_pdc) = c("Variables", "Partial Deviance Change")

EVI_pdc

EVI_dev = round(summary(EVI_full)$dev.expl, digits = 4)*100


#Create dataframe of gam full summary info
EVI_df = as.data.frame(summary(EVI_full)$s.table)
attributes(EVI_df)$row.names = "NULL"

#Join dataframes together and create a table 
EVI_table = cbind(EVI_pdc, EVI_df) %>%
  arrange(desc(`Partial Deviance Change`))%>%
  kable(caption = "EVI standardized anomaly GAM model, Deviance Explained = 10.6%")

EVI_table %>%
  kableExtra::kable_styling() %>%
  kableExtra::save_kable(file = "/home/kjamerson/pedotransfer/Research/Figures/EVI_anomaly_TABLE_MID.html")



