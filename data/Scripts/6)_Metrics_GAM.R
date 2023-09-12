#Correlation plots of NIRv by vegetation cover
#Date Modified: 3/20/2023
#Authors: Kayla Jamerson and Zach Hoylman

#This script generates GAM models for the response variables of NIRv, EVI, and SIF and the predictor variables of soil water potential, VWC, VPD and soil temperature.
#This script has code for full Full GAM models, partial effect plots, partial deviance change tables, predicted vs observed plots and aggregated by condition for independent plots.


Sys.setenv("OPENSSL_CONF"="/dev/null")




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
library(statar)
# esquisser for plot help

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
  drop_na() 

EVI_binded = list.files('/home/kjamerson/pedotransfer/Research/EVI', full.names = T) %>%
  purrr::map(., read_csv) %>%
  bind_rows() %>%
  drop_na() 


SIF = read_csv('/home/kjamerson/pedotransfer/Research/mesonet_sif.csv') %>%
  select(station, time, sif) %>%
  mutate(year = lubridate::year(time))
colnames(SIF) = c("station", "date", 'sif', 'year')


#Import soil data
soils_binded = list.files(path = '~/pedotransfer/Research/Water_Potential/Avg_Soil', full.names = T) %>%
  purrr::map(.,read_csv) %>%
  bind_rows() %>%
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


#Merge NIRv data and biophysical data

NIRv = NIRv_binded %>%
  select(station, name, time, NIRv, NIRv_rolling) %>%
  mutate( year = lubridate::year(time)) %>%
  filter(year %in% c(2017:2021))
colnames(NIRv) = c("station", "name", "date", 'NIRv',"NIRv_rolling",'year')


EVI = EVI_binded %>%
  select(station, name, date, value) %>%
  mutate( year = lubridate::year(date)) 
colnames(EVI) = c("station", "name", "date", 'EVI', 'year')



#Get Soils Data ready for merging. # Select mid soils
mid_Soils = soils_binded %>%  
  select(station, date, psi_kPa_mid, vpd_kPa_mid, vwc_mid, soilt_C_mid)%>%
  mutate(year = lubridate::year(date)) 


#Merge soil and NIR
NIRv_Soils = left_join(mid_Soils, NIRv, by = c('station' = 'station', 'date' = 'date', 'year'='year'))

EVI_Soils = left_join(mid_Soils, EVI, by = c('station' = 'station', 'date' = 'date', 'year'='year'))


SIF_Soils = left_join(mid_Soils, SIF, by = c('station' = 'station', 'date' = 'date', 'year'='year')) %>%
  drop_na()

#Merge with RAP Data
EVI_Data = left_join(EVI_Soils, RAP_Data, by = c('station' = 'station', 'year' = 'year')) %>%
  select(-c('name.x'))   

NIRv_Data = left_join(NIRv_Soils, RAP_Data, by = c('station' = 'station', 'year' = 'year')) %>%
  select(-c('name.x'))

SIF_Data = left_join(SIF_Soils, RAP_Data, by = c('station' = 'station', 'year' = 'year')) 



#Land_Cover = RAP_Data$top_groups %>%
 # unique()

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Filter out bad SWRC data
valid_sites = MSE_data %>%
  filter(MSE_data$Depth == '08',
         MSE < 0.00015)

#Summarize weekly data for all vegetation types, filtered for the growing season
EVI_Final = EVI_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(EVI > 0) %>%
  mutate(station_year = paste0(station,'_','year')) %>%
  filter(station_year %notin% mesonet$name_start_year,
      top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  mutate(month = lubridate::month(date),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         psi_kPa_mid_sqrt = sqrt(psi_kPa_mid))%>%
  filter(month %in% c(5:9)) %>%
  group_by(week, year, station) %>%
  summarise(psi_kPa_mid = median(psi_kPa_mid),
            psi_kPa_mid_sqrt = median(psi_kPa_mid_sqrt),
            psi_kPa_mid_log = median(psi_kPa_mid %>% log()),
            vpd_kPa_mid = median(vpd_kPa_mid),
            vwc_mid = median(vwc_mid),
            soilt_C_mid = median(soilt_C_mid),
            n = length(EVI),
            EVI = median(EVI))



NIRv_Final = NIRv_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(NIRv > 0) %>%
  mutate(station_year = paste0(station,'_','year')) %>%
  filter(station_year %notin% mesonet$name_start_year,
         top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  mutate(month = lubridate::month(date),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         psi_kPa_mid_sqrt = sqrt(psi_kPa_mid))%>%
  filter(month %in% c(5:9)) %>%
  group_by(week, year, station) %>%
  summarise(psi_kPa_mid = median(psi_kPa_mid),
            psi_kPa_mid_sqrt = median(psi_kPa_mid_sqrt),
            psi_kPa_mid_log = median(psi_kPa_mid %>% log()),
            vpd_kPa_mid = median(vpd_kPa_mid),
            vwc_mid = median(vwc_mid),
            soilt_C_mid = median(soilt_C_mid),
            n = length(NIRv),
            NIRv = median(NIRv)) 



SIF_Final = SIF_Data %>%
  filter(station %in% valid_sites$station,
         sif > 0) %>%
  mutate(station_year = paste0(station,'_','year')) %>%
  filter(station_year %notin% mesonet$name_start_year,
        top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  mutate(month = lubridate::month(date),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         psi_kPa_mid_sqrt = sqrt(psi_kPa_mid))%>%
  filter(month %in% c(5:9)) %>%
  group_by(week, year, station) %>%
  summarise(psi_kPa_mid = median(psi_kPa_mid),
            psi_kPa_mid_sqrt = median(psi_kPa_mid_sqrt),
            psi_kPa_mid_log = median(psi_kPa_mid %>% log()),
            vpd_kPa_mid = median(vpd_kPa_mid),
            vwc_mid = median(vwc_mid),
            soilt_C_mid = median(soilt_C_mid),
            n = length(sif),
            sif = median(sif)) 

# Fit full GAM model with all predictors 
k_ = 6

EVI_full <- gam(EVI ~ s(psi_kPa_mid_sqrt, k = k_) + s(vpd_kPa_mid, k = k_) + s(vwc_mid, k = k_) + s(soilt_C_mid, k = k_), family = Gamma(link = 'log'), data = EVI_Final, dispersion = 0.5)
NIRv_full <- gam(NIRv ~ s(psi_kPa_mid_sqrt, k = k_) + s(vpd_kPa_mid, k = k_) + s(vwc_mid, k = k_) + s(soilt_C_mid, k = k_), family = Gamma(link = 'log'), data = NIRv_Final, dispersion = 0.5)
SIF_full <- gam(sif ~ s(psi_kPa_mid_sqrt, k = k_) + s(vpd_kPa_mid, k = k_) + s(vwc_mid, k = k_) + s(soilt_C_mid, k = k_), family = Gamma(link = 'log'), data = SIF_Final, dispersion = 0.5)

EVI_rse = mean(summary(EVI_full$residuals^2))
NIRv_rse = mean(summary(NIRv_full$residuals^2))
SIF_rse = mean(summary(SIF_full$residuals^2))


#Create Partial Effect plots for each predictor variable 
pe_psi_NIRv = smooth_estimates(NIRv_full, smooth = "s(psi_kPa_mid_sqrt)") %>%
  mutate(Model = "NIRv")

pe_vpd_NIRv = smooth_estimates(NIRv_full, smooth = "s(vpd_kPa_mid)") %>%
  mutate(Model = "NIRv")

pe_vwc_NIRv = smooth_estimates(NIRv_full, smooth = "s(vwc_mid)") %>%
  mutate(Model = "NIRv")

pe_soilt_NIRv = smooth_estimates(NIRv_full, smooth = "s(soilt_C_mid)") %>%
  mutate(Model = "NIRv")

pe_psi_EVI = smooth_estimates(EVI_full, smooth = "s(psi_kPa_mid_sqrt)") %>%
  mutate(Model = "EVI")


pe_vpd_EVI = smooth_estimates(EVI_full, smooth = "s(vpd_kPa_mid)") %>%
  mutate(Model = "EVI")

pe_vwc_EVI = smooth_estimates(EVI_full, smooth = "s(vwc_mid)") %>%
  mutate(Model = "EVI")

pe_soilt_EVI = smooth_estimates(EVI_full, smooth = "s(soilt_C_mid)") %>%
  mutate(Model = "EVI")

pe_psi_SIF = smooth_estimates(SIF_full, smooth = "s(psi_kPa_mid_sqrt)") %>%
  mutate(Model = "SIF")

pe_vpd_SIF = smooth_estimates(SIF_full, smooth = "s(vpd_kPa_mid)") %>%
  mutate(Model = "SIF")

pe_vwc_SIF = smooth_estimates(SIF_full, smooth = "s(vwc_mid)") %>%
  mutate(Model = "SIF")

pe_soilt_SIF = smooth_estimates(SIF_full, smooth = "s(soilt_C_mid)") %>%
  mutate(Model = "SIF")

psi_pe_data = bind_rows(pe_psi_NIRv, pe_psi_EVI, pe_psi_SIF)
vpd_pe_data = bind_rows(pe_vpd_NIRv, pe_vpd_EVI, pe_vpd_SIF)
vwc_pe_data = bind_rows(pe_vwc_NIRv, pe_vwc_EVI, pe_vwc_SIF)
soilt_pe_data = bind_rows(pe_soilt_NIRv, pe_soilt_EVI, pe_soilt_SIF)


pe_psi_plot = psi_pe_data %>%
  add_confint() %>%
  ggplot(aes(y = est , x = psi_kPa_mid_sqrt, fill = Model))+
  theme_bw(base_size = 20)+
geom_line(aes(y = est , x = psi_kPa_mid_sqrt, color = Model), size = 1)+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.1) +
  labs(y = "Partial Effect",
       title = expression("Partial effect of predictor variables from GAM models"),
       x = expression("Ψ (- sqrt kPa)"), element_text( size = 25))

pe_psi_plot

pe_vpd_plot = vpd_pe_data %>%
  add_confint() %>%
  ggplot(aes(y = est , x = vpd_kPa_mid, fill = Model))+
  theme_bw(base_size = 20)+
  geom_line(aes(y = est , x = vpd_kPa_mid, color = Model), size = 1)+
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.1) +
  labs(y = "Partial Effect",
       x = expression("VPD (kPa)"), element_text( size = 25))

pe_vpd_plot

pe_vwc_plot = vwc_pe_data %>%
  add_confint() %>%
  ggplot(aes(y = est , x = vwc_mid, fill = Model))+
  theme_bw(base_size = 20)+
  geom_line(aes(y = est , x = vwc_mid, color = Model), size = 1)+
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.1) +
  labs(y = "Partial Effect",
       x = expression("VWC (m³/m³)"), element_text( size = 25))

pe_vwc_plot

pe_soilt_plot = soilt_pe_data %>%
  add_confint() %>%
  ggplot(aes(y = est , x = soilt_C_mid, fill = Model))+
  theme_bw(base_size = 20)+
  geom_line(aes(y = est , x = soilt_C_mid, color = Model), size = 1)+
   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.1) +
  labs(y = "Partial Effect",
       x = expression("Soil Temp (C)"), element_text( size = 25))

pe_soilt_plot

PE_Plots_All = (pe_psi_plot / pe_vpd_plot) | (pe_vwc_plot / pe_soilt_plot)
PE_Plots_All
ggsave(PE_Plots_All, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/ALL_PE_Plots_in.png'), width = 28, height = 14, units = 'in')



#Fit a null model without each predictor variable for calculating partial deviance change  
NIRv_null_psi_kPa_mid = gam(NIRv ~  s(vpd_kPa_mid) + s(vwc_mid) + s(soilt_C_mid), sp = NIRv_full$sp[-1], family = Gamma, data = NIRv_Final, dispersion = 0.5)
NIRv_null_vpd_kPa_mid = gam(NIRv ~ s(psi_kPa_mid_sqrt) +  s(vwc_mid) + s(soilt_C_mid), sp = NIRv_full$sp[-2], family = Gamma, data = NIRv_Final, dispersion = 0.5)
NIRv_null_vwc_mid = gam(NIRv ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid)  + s(soilt_C_mid), sp = NIRv_full$sp[-3], family = Gamma, data = NIRv_Final, dispersion = 0.5)
NIRv_null_soilt_C_mid = gam(NIRv ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid) + s(vwc_mid), sp = NIRv_full$sp[-4], family = Gamma, data = NIRv_Final, dispersion = 0.5)

EVI_null_psi_kPa_mid = gam(EVI ~  s(vpd_kPa_mid) + s(vwc_mid) + s(soilt_C_mid), sp = EVI_full$sp[-1], family = Gamma, data = EVI_Final, dispersion = 0.5)
EVI_null_vpd_kPa_mid = gam(EVI ~ s(psi_kPa_mid_sqrt) +  s(vwc_mid) + s(soilt_C_mid), sp = EVI_full$sp[-2], family = Gamma, data = EVI_Final, dispersion = 0.5)
EVI_null_vwc_mid = gam(EVI ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid)  + s(soilt_C_mid), sp = EVI_full$sp[-3], family = Gamma, data = EVI_Final, dispersion = 0.5)
EVI_null_soilt_C_mid = gam(EVI ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid) + s(vwc_mid), sp = EVI_full$sp[-4], family = Gamma, data = EVI_Final, dispersion = 0.5)

SIF_null_psi_kPa_mid = gam(sif ~  s(vpd_kPa_mid) + s(vwc_mid) + s(soilt_C_mid), sp = SIF_full$sp[-1], family = Gamma, data = SIF_Final, dispersion = 0.5)
SIF_null_vpd_kPa_mid = gam(sif ~ s(psi_kPa_mid_sqrt) +  s(vwc_mid) + s(soilt_C_mid), sp = SIF_full$sp[-2], family = Gamma, data = SIF_Final, dispersion = 0.5)
SIF_null_vwc_mid = gam(sif ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid)  + s(soilt_C_mid), sp = SIF_full$sp[-3], family = Gamma, data = SIF_Final, dispersion = 0.5)
SIF_null_soilt_C_mid = gam(sif ~ s(psi_kPa_mid_sqrt) + s(vpd_kPa_mid) + s(vwc_mid), sp = SIF_full$sp[-4], family = Gamma, data = SIF_Final, dispersion = 0.5)

# Subtract the null models from the full models. 
#NIRv
a_NIRv = round((deviance(NIRv_null_psi_kPa_mid)-deviance(NIRv_full)), digits = 5)
b_NIRv = round((deviance(NIRv_null_vpd_kPa_mid)-deviance(NIRv_full)), digits =5)
c_NIRv = round((deviance(NIRv_null_vwc_mid)-deviance(NIRv_full)), digits = 5)
d_NIRv = round((deviance(NIRv_null_soilt_C_mid)-deviance(NIRv_full)), digits = 5)

#EVI
a_EVI = round((deviance(EVI_null_psi_kPa_mid)-deviance(EVI_full)), digits = 5)
b_EVI = round((deviance(EVI_null_vpd_kPa_mid)-deviance(EVI_full)), digits =5)
c_EVI = round((deviance(EVI_null_vwc_mid)-deviance(EVI_full)), digits = 5)
d_EVI = round((deviance(EVI_null_soilt_C_mid)-deviance(EVI_full)), digits = 5)

#SIF
a_SIF = round((deviance(SIF_null_psi_kPa_mid)-deviance(SIF_full)), digits = 5)
b_SIF = round((deviance(SIF_null_vpd_kPa_mid)-deviance(SIF_full)), digits =5)
c_SIF= round((deviance(SIF_null_vwc_mid)-deviance(SIF_full)), digits = 5)
d_SIF = round((deviance(SIF_null_soilt_C_mid)-deviance(SIF_full)), digits = 5)

#Create dataframe of partial Deviance Change
predictor_names <- c( "Ψ (sqrt kPa)","VPD (kPa)","VWC (m³/m³)", "Soil Temperature (C)")

pdc_NIRv = data.frame( variables = predictor_names, pdc_NIRv = c(a_NIRv,b_NIRv,c_NIRv,d_NIRv))
colnames(pdc_NIRv) = c("Variables", "Partial Deviance Change")
pdc_NIRv %>% arrange(desc(`Partial Deviance Change`))

pdc_EVI = data.frame( variables = predictor_names, pdc_EVI = c(a_EVI,b_EVI,c_EVI,d_EVI))
colnames(pdc_EVI) = c("Variables", "Partial Deviance Change")
pdc_EVI %>% arrange(desc(`Partial Deviance Change`)) 

pdc_SIF = data.frame(variables = predictor_names, pdc_SIF = c(a_SIF,b_SIF,c_SIF,d_SIF))
colnames(pdc_SIF) = c("Variables", "Partial Deviance Change")
pdc_SIF %>% arrange(desc(`Partial Deviance Change`))

NIRv_dev = round(summary(NIRv_full)$dev.expl, digits = 4)*100

#Create dataframe of gam full summary info
NIRv_df = as.data.frame(summary(NIRv_full)$s.table)
attributes(NIRv_df)$row.names = "NULL"

#Join dataframes together and create a table 
NIRv_table = cbind(pdc_NIRv, NIRv_df) %>%
  mutate(Model = "NIRv", .before = "Variables") %>%
  mutate(Deviance_Explained = "20.14%") %>%
  arrange(desc(`Partial Deviance Change`)) 
colnames(NIRv_table)[8]= "Deviance Explained"



EVI_dev = round(summary(EVI_full)$dev.expl, digits = 4)*100

#Create dataframe of gam full summary info
EVI_df = as.data.frame(summary(EVI_full)$s.table)
attributes(EVI_df)$row.names = "NULL"

#Join dataframes together and create a table 

EVI_table = cbind(pdc_EVI, EVI_df) %>%
  mutate(Model = "EVI", .before = "Variables") %>%
  mutate(Deviance_Explained = "24.54%") %>%
  arrange(desc(`Partial Deviance Change`)) 
colnames(EVI_table)[8]= "Deviance Explained"


SIF_dev = round(summary(SIF_full)$dev.expl, digits = 4)*100

SIF_df = as.data.frame(summary(SIF_full)$s.table)
attributes(SIF_df)$row.names = "NULL"

#Join dataframes together and create a table 

SIF_table = cbind(pdc_SIF, SIF_df) %>%
  mutate(Model = "SIF", .before = "Variables") %>%
  mutate(Deviance_Explained = "24.93%") %>%
  arrange(desc(`Partial Deviance Change`)) 
colnames(SIF_table)[8]= "Deviance Explained"


All_tables = bind_rows(NIRv_table, EVI_table, SIF_table) %>%
  kable() %>%
  kableExtra::kable_styling() %>%
  kableExtra::row_spec(1, 1, extra_css = "border-bottom: 1px solid")%>%
  kableExtra::row_spec(5, 5, extra_css = "border-bottom: 1px solid")%>%
  kableExtra::row_spec(9, 9, extra_css = "border-bottom: 1px solid")%>%
  kableExtra::save_kable(file = '//home/kjamerson/pedotransfer/Research/Final_Figures/ALL_table.png')

#predicted vs observed plots
NIRv_psi_kPa_pred = 
  ggplot()+
  labs(x = 'Predicted', y = 'Observed')+
  geom_point(data = NULL, aes(x = predict(NIRv_full, type = 'response'), y = NIRv_Final$NIRv),size = 3)+
  scale_color_discrete(name = NULL)+
  theme_bw(base_size = 14)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_smooth(method = 'lm')+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  ggtitle("NIRv Predicted vs. Observed")+
  theme(plot.title = element_text( size = 27),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))

NIRv_psi_kPa_pred

EVI_psi_kPa_pred = 
  ggplot()+
  labs(x = 'Predicted', y = 'Observed')+
  geom_point(data = NULL, aes(x = predict(EVI_full, type = 'response'), y = EVI_Final$EVI),size = 3)+
  scale_color_discrete(name = NULL)+
  theme_bw(base_size = 14)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_smooth(method = 'lm')+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  ggtitle("EVI Predicted vs. Observed")+
  theme(plot.title = element_text( size = 27),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))

EVI_psi_kPa_pred

SIF_psi_kPa_pred = 
  ggplot()+
  labs(x = 'Predicted', y = 'Observed')+
  geom_point(data = NULL, aes(x = predict(SIF_full, type = 'response'), y = SIF_Final$sif),size = 3)+
  scale_color_discrete(name = NULL)+
  theme_bw(base_size = 14)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  geom_smooth(method = 'lm')+
  ggtitle("SIF Predicted vs. Observed")+
  geom_abline(intercept = 0, slope = 1, color = 'red')+
  theme(plot.title = element_text( size = 27),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))

SIF_psi_kPa_pred

Predicted_Plots = cowplot::plot_grid(NIRv_psi_kPa_pred, EVI_psi_kPa_pred, SIF_psi_kPa_pred, nrow = 3,  labels = c('A', 'B', 'C'))
ggsave(Predicted_Plots, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/Predicted_Plots_in.png'), width = 9, height = 14, units = 'in')

# Aggregate by condition
n_out = 50

psi_aggregation = NIRv_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(NIRv > 0) %>%
  mutate(station_year = paste0(station,'_','year'),
         month = lubridate::month(date)) %>%
  filter(station_year %notin% mesonet$name_start_year,#) %>%#,
         top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  filter(month %in% c(5:9)) %>%
  mutate(bin = .bincode(psi_kPa_mid %>% sqrt, breaks = seq(0.1, max(psi_kPa_mid) %>% sqrt, length.out = n_out))) %>%
  group_by(bin) %>%
  summarise(NIRv_median = median(NIRv),
            NIRv_upper = quantile(NIRv, 0.75),
            NIRv_lower = quantile(NIRv, 0.25),
            NIRv_IQR = NIRv_upper - NIRv_lower,
            psi = median(psi_kPa_mid %>% sqrt))

psi_gam = gam(NIRv_median ~ s(psi, k = k_), weights = 1-NIRv_IQR, data = psi_aggregation, family = Gamma(link = 'log')) 

psi_fit = tibble(predicted_NIRv = 
                   predict(psi_gam, type = 'response', newdata = 
                             data.frame(psi = seq(min(psi_aggregation$psi), max(psi_aggregation$psi), length.out = n_out))),
                 psi = seq(min(psi_aggregation$psi), max(psi_aggregation$psi), length.out = n_out)
)

psi_ag_plot = psi_aggregation %>%
  ggplot() + 
  geom_ribbon(aes(x = psi, y = NIRv_median, ymax = NIRv_upper, ymin = NIRv_lower), alpha = 0.1)+
  geom_point(aes(x = psi, y = NIRv_median), size = 2)+
  labs(x = 'Ψ (- sqrt kPa)', y = 'Median NIRv')+
  ggtitle("Aggregation by condition ")+
  geom_line(data = psi_fit, aes(x = psi, y = predicted_NIRv), color = "blue") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text( size = 27),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))

#VPD 
vpd_aggregation = NIRv_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(NIRv > 0) %>%
  mutate(station_year = paste0(station,'_','year'),
         month = lubridate::month(date)) %>%
  filter(station_year %notin% mesonet$name_start_year,
         top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  filter(month %in% c(5:9)) %>%
  mutate(bin = .bincode(vpd_kPa_mid , breaks = seq(0.1, quantile(vpd_kPa_mid, 0.95), length.out = n_out))) %>%
  group_by(bin) %>%
  summarise(NIRv_median = median(NIRv),
            NIRv_upper = quantile(NIRv, 0.75),
            NIRv_lower = quantile(NIRv, 0.25),
            NIRv_IQR = NIRv_upper - NIRv_lower,
            vpd = median(vpd_kPa_mid))

vpd_gam = gam(NIRv_median ~ s(vpd, k = k_), weights = 1-NIRv_IQR, data = vpd_aggregation, family = Gamma(link = 'log')) 

vpd_fit = tibble(predicted_NIRv = 
                   predict(vpd_gam, type = 'response', newdata = 
                             data.frame(vpd = seq(min(vpd_aggregation$vpd), max(vpd_aggregation$vpd), length.out = n_out))),
                 vpd = seq(min(vpd_aggregation$vpd), max(vpd_aggregation$vpd), length.out = n_out)
)

vpd_ag_plot = vpd_aggregation %>%
  ggplot() + 
  geom_ribbon(aes(x = vpd, y = NIRv_median, ymax = NIRv_upper, ymin = NIRv_lower), alpha = 0.1)+
  geom_point(aes(x = vpd, y = NIRv_median), size = 2)+
  labs(x = 'VPD (kPa)', y = 'Median NIRv')+
  theme_bw(base_size = 14) +
  geom_line(data = vpd_fit, aes(x = vpd, y = predicted_NIRv), color = "blue")+
  theme(plot.title = element_text( size = 27),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))

# VWC
vwc_aggregation = NIRv_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(NIRv > 0) %>%
  mutate(station_year = paste0(station,'_','year'),
         month = lubridate::month(date)) %>%
  filter(station_year %notin% mesonet$name_start_year,
         top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  filter(month %in% c(5:9)) %>%
  mutate(bin = .bincode(vwc_mid , breaks = seq(0.1, quantile(vwc_mid, 0.95), length.out = n_out))) %>%
  group_by(bin) %>%
  summarise(NIRv_median = median(NIRv),
            NIRv_upper = quantile(NIRv, 0.75),
            NIRv_lower = quantile(NIRv, 0.25),
            NIRv_IQR = NIRv_upper - NIRv_lower,
            vwc = median(vwc_mid))

vwc_gam = gam(NIRv_median ~ s(vwc, k = k_), weights = 1-NIRv_IQR, data = vwc_aggregation, family = Gamma(link = 'log')) 

vwc_fit = tibble(predicted_NIRv = 
                   predict(vwc_gam, type = 'response', newdata = 
                             data.frame(vwc = seq(min(vwc_aggregation$vwc), max(vwc_aggregation$vwc), length.out = n_out))),
                 vwc = seq(min(vwc_aggregation$vwc), max(vwc_aggregation$vwc), length.out = n_out)
)

vwc_ag_plot = vwc_aggregation %>%
  ggplot() + 
  geom_ribbon(aes(x = vwc, y = NIRv_median, ymax = NIRv_upper, ymin = NIRv_lower), alpha = 0.1)+
  geom_point(aes(x = vwc, y = NIRv_median), size = 2)+
  labs(x = 'VWC (m³/m³)', y = 'Median NIRv')+
  theme_bw(base_size = 14) +
  geom_line(data = vwc_fit, aes(x = vwc, y = predicted_NIRv), color = "blue")+
  theme(plot.title = element_text( size = 27),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))



# Soil temp
soilt_aggregation = NIRv_Data %>%
  filter(station %in% valid_sites$station) %>%
  filter(NIRv > 0) %>%
  mutate(station_year = paste0(station,'_','year'),
         month = lubridate::month(date)) %>%
  filter(station_year %notin% mesonet$name_start_year,
         top_groups %notin% c(NA, "No Dominant Class", 'Other')) %>%
  filter(month %in% c(5:9)) %>%
  mutate(bin = .bincode(soilt_C_mid , breaks = seq(0.1, quantile(soilt_C_mid, 0.95), length.out = n_out))) %>%
  group_by(bin) %>%
  summarise(NIRv_median = median(NIRv),
            NIRv_upper = quantile(NIRv, 0.75),
            NIRv_lower = quantile(NIRv, 0.25),
            NIRv_IQR = NIRv_upper - NIRv_lower,
            soilt = median(soilt_C_mid))

soilt_gam = gam(NIRv_median ~ s(soilt, k = k_), weights = 1-NIRv_IQR, data = soilt_aggregation, family = Gamma(link = 'log')) 

soilt_fit = tibble(predicted_NIRv = 
                   predict(soilt_gam, type = 'response', newdata = 
                             data.frame(soilt = seq(min(soilt_aggregation$soilt), max(soilt_aggregation$soilt), length.out = n_out))),
                 soilt = seq(min(soilt_aggregation$soilt), max(soilt_aggregation$soilt), length.out = n_out)
)

soilt_ag_plot = soilt_aggregation %>%
  ggplot() + 
  geom_ribbon(aes(x = soilt, y = NIRv_median, ymax = NIRv_upper, ymin = NIRv_lower), alpha = 0.1)+
  geom_point(aes(x = soilt, y = NIRv_median), size = 2)+
  labs(x = 'Soil Temp (C)', y = 'Median NIRv')+
  theme_bw(base_size = 14) +
  geom_line(data = soilt_fit, aes(x = soilt, y = predicted_NIRv), color = "blue")+
  theme(plot.title = element_text( size = 27),
        legend.text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25))


Aggregation_Plots = cowplot::plot_grid(psi_ag_plot, vpd_ag_plot, vwc_ag_plot, soilt_ag_plot, nrow = 2)
ggsave(Aggregation_Plots, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/Aggregation_Plots_in.png'), width = 14, height = 14, units = 'in')



