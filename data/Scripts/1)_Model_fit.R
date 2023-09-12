
#Water Retention Characterization
#Date Modified: 3/13/2023
#Authors: Zach Hoylman, modified by Kayla Jamerson

#This script is designed for data that has been exported from HYPROP Fit software following HYPROP laboratory analysis. It generates soil water retention models using
# the Fredlund-Xing Model,Van Genechten Model,Kosugi Model, and the Brooks Model. The Mean Standard Error is then used to select the best model with the best fit. 

#Install packages


#Load libraries 
library(ggplot2)
require(scales)
library(readxl)   
library(tcltk)
library(grid)
library(gridExtra)
library(cowplot)
library(minpack.lm)
library(dplyr)
library(readr)
library(stringr)
library(HydroMe)
library(magrittr)
library(dplyr)
library(tidyr)
library(pracma)
library(SciViews)


#Create empty dataframe to store MSE values in 
x = c("station", "Depth", "MSE", "model")
MSE_final = data.frame((matrix(nrow = 0, ncol = length(x))))
colnames(MSE_final) = x


#Create function to read in data 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Create logarithmic spaced sequence 
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#Write soils function to fit data 
fit_soils = function(data){
  file_name = basename(data)
  data = read_excel_allsheets(data)
  
  #select sheet
  pf_vwc = data$`Evaluation-Retention Θ(pF)`%>%
    dplyr::filter(`pF [-]` >= 1)
  
  
  #Convert pF into Pa
  MPa = ((10^(pf_vwc$`pF [-]`))/10200)
  kPa = MPa * 1000
  
  #Convert theta into ratio
  VWC = pf_vwc$`Water Content [Vol%]`/100
  
  #Fit models
 # van_model <- nlsLM(VWC ~ r + (s - r)/((1+a*(abs(kPa)^n))^(1-(1/n))), start = list(r = 0, s = 0.483, a = 0.1966, n = 1.241))
  van_model <- nlsLM(VWC ~ r + (s - r)*((1/(1+a*(abs(kPa)^n)))^(1-(1/n))), start = list(r = 0, s = 0.483, a = 0.1966, n = 1.241))


  brooks_model <- 
    tryCatch({
    nlsLM(VWC ~ r + ((s-r)*((kPa/a)^(-n))),
                       start = list(r = 0, s = 0.483, a = 0.1966, n = 1.241), 
                      control = nls.lm.control(maxiter=500))      
  }, error = function(e){
   brooks_model <- nlsLM(VWC ~ r + ((s-r)*((kPa/a)^(-n))),
          start = list(r = 0, s = 0.0001, a = 0.0001, n = 0.0001), 
          control = nls.lm.control(maxiter=500))    
   return(brooks_model) 
  })
  
  Kosugi_model<- 
   tryCatch({
    nlsLM(VWC ~ r + (0.5)*(s-r)*erfc(((ln(kPa/m))/(sigma*sqrt(2)))), start = list(r = 0, s = 0.483, sigma = 0.1966, m = 1.241),
          control = nls.lm.control(maxiter=500))
 
   }, warning = function(w){
    Kosugi_model <- nlsLM(VWC ~ r + (0.5)*(s-r)*erfc(((ln(kPa/m))/(sigma*sqrt(2)))), start = list(r = 0, s = 0.483, sigma = 0.1966, m = 1.241),
           control = nls.lm.control(maxiter=500))
   
    }, error = function(e){
      Kosugi_model<-  nlsLM(VWC ~ r + (0.5)*(s-r)*erfc(((ln(kPa/m))/(sigma*sqrt(2)))), start = list(r = 0, s = 0.01, sigma = 0.1966, m = 1.241),
                            control = nls.lm.control(maxiter=500))
      return(Kosugi_model)
    }, error = function(e){
      Kosugi_model<-  nlsLM(VWC ~ r + (0.5)*(s-r)*erfc(((ln(kPa/m))/(sigma*sqrt(2)))), start = list(r = 0, s = 0.01, sigma = 0.8, m = 1.241),
                            control = nls.lm.control(maxiter=500))
      return(Kosugi_model) 
    })



  FX_model <-
    tryCatch ({
    nlsLM(VWC ~ r + (s-r)/((log(exp(1)+((kPa/h)^(n)))^(m))), start = list( r = 0, s = 0.483, n = 1.241, m=0.559, h =0.01),
          control = nls.lm.control(maxiter=500))
      }, error = function(e){
  FX_model = nlsLM(VWC ~ r + (s-r)/((log(exp(1)+((kPa/h)^(n)))^(m))), start = list( r = 0.8, s = 0.483, n = 1.241, m=0.559, h =0.01))
  return(FX_model)
})
  


  #Store model names
  models = list(FX_model, van_model, Kosugi_model, brooks_model)
  model_names = c("Fredlund-Xing Model","Van Genechten Model","Kosugi Model","Brooks Model")
  
  FX_summary <- summary(FX_model)
  
  van_summary <- summary(van_model)
 
 Kosugi_summary <-summary(Kosugi_model)
 
  brooks_summary <-  summary(brooks_model)

  #Select Best Fit 
  FX_MSE <- mean(FX_summary$residuals^2)
 van_MSE <- mean(van_summary$residuals^2)
Kosugi_MSE <- mean(Kosugi_summary$residuals^2) 
 brooks_MSE <- mean(brooks_summary$residuals^2)
  
  MSE_df <- data_frame(FX_MSE, van_MSE, Kosugi_MSE, brooks_MSE, na.rm = TRUE)
  
  Best_Fit <- case_when(
    min(MSE_df) == MSE_df[1] ~ "Fredlund-Xing Model",
    min(MSE_df) == MSE_df[2] ~ "Van Genechten Model",
    min(MSE_df) == MSE_df[3] ~ "Kosugi Model",
    min(MSE_df) == MSE_df [4] ~ "Brooks Model")
  
  best_model = models[[which(model_names == Best_Fit)]] 
  
  dummy_kPa = data.frame(kPa = lseq(min(kPa), max(kPa), 10000))
  
  #Predict theoretical curve using dummy data
  model_predict = predict(best_model, newdata = dummy_kPa)

  #Store model coeffitients for inverse function
  model_coef = coef(best_model)

  
  #compile export list
  export = list()
  export[[1]] = best_model
  export[[2]] = data.frame(raw_kPa = kPa,
                           raw_VWC = VWC)
  export[[3]] = data.frame(kPa = dummy_kPa,
                           fit_VWC = model_predict)
  export[[4]] = file_name
  export[[5]] = model_names[which(model_names == Best_Fit)]
  return(export)
  
}


#Create water rentetion plot
plot_data = function(data){
  
  stations = substr(model[4], 1, 8) 
  
  depth = substr(model[4], 9,10)
  
  depth_title <- case_when(
    depth == "02" ~ "4 Inch Depth",
    depth == "04" ~ "4 Inch Depth",
    depth == "08" ~ "8 Inch Depth",
    depth == "20" ~ "20 Inch Depth")
  
 
  
  ggplot_data = data.frame(kPa = data[[2]]$raw_kPa, VWC = data[[2]]$raw_VWC)
  
  model_data = data.frame(model_predict = data[[3]]$fit_VWC, dummy_kPa = data[[3]]$kPa)
  
 SWRC_Plot = ggplot(data = ggplot_data, aes(x =VWC , y = kPa))+
    geom_point(size = 3, color = "#BFEFFF")+
    geom_point(size = 3, color = "black", fill = NA, shape = 21, alpha = 0.3)+
    theme_bw()+
    scale_y_log10(labels = comma)+
    scale_x_continuous(breaks = pretty(ggplot_data$VWC, n = 10))+
    ylab("Ψ (kPa)")+
    xlab(expression(paste("Volumetric Water Content (", m^3," ", m^-3,")")))+
    geom_line(data = model_data, aes(x = model_predict, y = dummy_kPa))+
    theme(text = element_text(size=18),
          plot.margin=unit(c(2,2,2,2),"cm"),
          plot.caption = element_text(hjust = 0.5),
          axis.line = element_line(color='black'),
          panel.grid.major = element_line(colour="grey", size=0.25))+
    geom_hline(yintercept=c(1500), linetype = "dashed", color = "red")+
    annotate(geom = "text", x = (max(ggplot_data$VWC)*.7), y = 3000, hjust = 1.05,
             label = "Theoretical Wilting Point", size = 5, color = "red")+
    labs(title = stations,
            subtitle = depth_title, 
            caption = data[[5]]) 
  

  SWRC_Plot
  
}

#Select folder where HYPROP data is stored 
work.dir = '/home/kjamerson/pedotransfer/Research/Model_Data'

data = list.files(work.dir, pattern = ".xlsx$", full.names = T)


# select write dir for plot export
write.dir = '/home/kjamerson/pedotransfer/Research/SWRC_Plots'

#Out_data information to view
out_data = data.frame(mse = rep(NA, length(data)),
                      model = rep(NA, length(data)),
                      depth = substr(data, 61, 62))

for(i in 1:length(data)){


  tictoc::tic()
  #run model fit
  model = fit_soils(data[i])
  
  #plot and export 
  plot_data(model)

 
   ggsave(filename = paste0(write.dir,'/',model[[4]],"_water_retension_curve.png"), plot = last_plot(), dpi = 500,
         units = "in", width = 10, height = 6)
  
  saveRDS(model, paste0(write.dir,'/', model[[4]], "_model.RDS"))
  capture.output(summary(model[[1]]), file = paste0(write.dir,'/', model[[4]], "_model_information.txt"))
  tictoc::toc()
  
  stations = substr(model[4], 1, 8) 
  
  depth = substr(model[4], 9,10)
  
  out_data$model[i] = model[[5]]
  
  out_data$mse[i] = mean(summary(model[[1]])$residuals^2)
    


  #MSE dataframe  
MSE = data.frame(
  "station" = stations,
  "Depth" = depth,
  "MSE" = out_data$mse[i],
  "model" = out_data$model[i])



#Bind all MSE data together   
MSE_final = rbind(MSE_final, MSE)
  
  

  print(i)
  
  
}

#Save out MSE Data for binding in Script 4
models = MSE_final %>%
  filter(Depth == "08") 
write_csv(MSE_final, paste0('/home/kjamerson/pedotransfer/Research/MSE_Model_Data'))


#View distribution of MSE data 
final = as.data.frame(models)
colnames(final) = c("Station","Depth", "MSE","Model")

table(final$model)

mse <- mse[order(as.numeric(out_data$mse))]

as.data.frame(MSE_final)

MSE_distribution =
 final %>%
  ggplot(mapping = aes(x = MSE, fill = Model)) +
  geom_histogram(color = "black")+
  geom_vline(xintercept = 0.00015, linetype = 2, color = "red")+
  labs (x = "Distribution of MSE values", y = "Count")+ #Create label
  theme_bw(base_size = 14) +
  theme(text=element_text(size = 18))+ #Make font size bigger
  scale_fill_discrete(labels = c("Fredlund-Xing (37)", "Kosugi (2)", "van Genuchten (13)"))+
  annotate(geom = "text", x =0.00024, y = 20,
           label = "0.00015", size = 8, color = "red")


ggsave(MSE_distribution, file = paste0('/home/kjamerson/pedotransfer/Research/Final_Figures/MSE_Distribution_in.png'), width = 10, height = 5, units = 'in')










