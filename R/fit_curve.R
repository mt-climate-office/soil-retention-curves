# Water Retention Characterization (Authors: Zach Hoylman, Kayla Jamerson, and Fin Malone)
require(scales)
library(readxl)
library(grid)
library(gridExtra)
library(cowplot)
library(minpack.lm)
library(pracma)
library(tidyverse)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, .name_repair = "unique_quiet"))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}


fit_soils = function(data){
  file_name = basename(data)
  data = read_excel_allsheets(data)
  
  # Select Sheet
  pF_vwc = data$`Evaluation-Retention Θ(pF)` %>% filter((`pF [-]` > 0) & (`Water Content [Vol%]` > 0))
  tryCatch({
    WP4C = data$`Evaluation-WP4 manual` %>%
      filter(!(`pF [-]` == 'NaN')) %>%
      mutate(`Water Content [Vol%]` = 
               ((as.numeric(`Gross wet mass [g]`) - as.numeric(`Tare mass [g]`)) -
                   (as.numeric(`Gross dry mass [g]`) - as.numeric(`Tare mass [g]`))) /
               (as.numeric(`Gross dry mass [g]`) - as.numeric(`Tare mass [g]`)) * 100 *
               data$Information$Value[which(data$Information$`Parameter Name` == 'Density [g/cm3]:')] %>% as.numeric(),
             
             `pF [-]` = as.numeric(`pF [-]`))
    
    pF_vwc = full_join(
      pF_vwc,
      WP4C %>% 
        select(`pF [-]`,`Water Content [Vol%]`),
      by=c("pF [-]", "Water Content [Vol%]")
    )
  })
  # Convert pF into kPa
  MPa = (10^(pF_vwc$`pF [-]`)) / 10000
  kPa = MPa * 1000
  
  # Convert Volumetric Water Content to fraction (m^-3/m^-3)
  VWC = pF_vwc$`Water Content [Vol%]`/100
  
  
  ######################## FIT MODELS
  
  # van Genuchten Model
  van <- nlsLM(VWC ~ O_r + (O_s - O_r)/((1+alpha*(abs(kPa)^n))^(1-(1/n))),
               start = list(O_r = 0, O_s = 0.483, alpha = 0.1966, n = 1.241),
               control = nls.lm.control(maxiter=500))
  
  
  # Kosugi Model
  Kosugi <- nlsLM( VWC ~ O_r + (O_s - O_r) * ((0.5*erfc(log(kPa/h_m) / (sqrt(2)*sigma)))),
                   control = nls.lm.control(maxiter=500),
                   start = list(O_r= 0, O_s = 0.483, sigma = 0.728, h_m = 116))
  
  
  # Fredlund-Xing Model
  
  FX = tryCatch({FX = nlsLM(VWC ~ O_r + (O_s - O_r)/((log(exp(1)+((kPa/h)^(n)))^(m))),           # Catches fitting error and changes
                            start = list( O_r = 0, O_s = 0.483, n = 1.241, m=0.559, h =0.01),    # O_r starting parameter
                            control = nls.lm.control(maxiter=500))},
                error = function(e)
                {FX = nlsLM(VWC ~ O_r + (O_s - O_r)/((log(exp(1)+((kPa/h)^(n)))^(m))),
                            start = list( O_r = 0.8, O_s = 0.483, n = 1.241, m=0.559, h =0.01),
                            control = nls.lm.control(maxiter=500))
                return(FX)})
  
  
  # Store Model Names
  models = list(van, Kosugi, FX) 
  model_names = c("Van Genuchten Model",
                  "Kosugi Model",
                  "Fredlund-Xing Model")
  van_summary <- summary(van)
  Kosugi_summary <- summary(Kosugi)
  FX_summary <- summary(FX)
  
  # Select Best Fit
  van_RSE <- van_summary$sigma
  Kosugi_RSE <- Kosugi_summary$sigma
  FX_RSE <- FX_summary$sigma
  
  
  RSE_df <- data_frame(van_RSE, Kosugi_RSE, FX_RSE)
  
  # Always use FX model because we have the inverse).  
  Best_Fit <- "Fredlund-Xing Model"
  best_model = models[[which(model_names == Best_Fit)]]
  
  
  dummy_kPa1 = data.frame(kPa = lseq(min(kPa), max(kPa), 10000))
  dummy_kPa2 = data.frame(kPa = lseq(.1, 6309573, 10000))
  
  # Predict Theoretical Curve Using Dummy Data
  model_predict1 = predict(best_model, newdata = dummy_kPa1)
  model_predict2 = predict(best_model, newdata = dummy_kPa2)
  
  # Store Model Coefficients for Inverse Function
  model_coef = coef(best_model)
  
  # Identify site name from file
  site_name = file_name
  
  # Compile Export List
  export = list()
  export[[1]] = best_model
  export[[2]] = data.frame(raw_kPa = kPa,
                           raw_VWC = VWC)
  export[[3]] = data.frame(kPa = dummy_kPa1,
                           fit_VWC = model_predict1)
  export[[4]] = file_name
  export[[5]] = model_names[which(model_names == Best_Fit)]
  export[[6]] = site_name
  export[[7]] = RSE_df
  export[[8]] = data.frame(kPa = dummy_kPa2,
                           fit_VWC = model_predict2)
  
  names(export) <- c("model", "raw_data", "fit_data", "filename", "best_fit", "filename2", "error", "fit_data2")
  return(export)
  
}

fit_all_depths <- function(zipped_data) {
  
  tmp <- tempdir()
  unzip(zipped_data, exdir=tmp)
  
  station <- basename(zipped_data) %>% 
    stringr::str_split_1("_") %>% 
    magrittr::extract(1)
  
  out <- list.files(
    file.path(tmp), pattern = '*.xlsx$', full.names = T, recursive=T
  ) %>% 
    purrr::map(function(x) {
      
      model = suppressWarnings(fit_soils(x))
      depth = stringr::str_extract(x, "\\d{2,3}cm")
      extra = stringr::str_detect(x, "_B.xlsx|_B_VC.xlsx")
      is_vc = stringr::str_detect(x, "_VC.xlsx")
      
      if (is.na(depth)) {
        # depth = str_extract(x, "(?<=_)(02|04|08|20|36)(?=_)") %>%
        depth = str_extract(x, "(?<=_)(02|04|08|20|36)(?=_)") %>%
          dplyr::case_match(
            # "02" ~ "10cm",
            "04" ~ "10cm",
            "08" ~ "20cm",
            "20" ~ "50cm", 
            "36" ~ "91cm"
          )
      }
      
      tibble::tibble(
        station = station,
        depth = depth,
        extra = extra,
        is_vc = is_vc,
        model = list(model), 
      )
    }) %>%
    dplyr::bind_rows() %>% 
    dplyr::group_by(depth) %>% 
    dplyr::filter(if (dplyr::n() > 1) (extra & is_vc) | is_vc else TRUE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-extra, -is_vc)
  
  unlink(file.path(tmp, "Post-Processing"), recursive = T)
  
  return(out)
}






