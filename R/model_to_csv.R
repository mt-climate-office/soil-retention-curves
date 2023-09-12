library(magrittr)

params_to_csv <- function(x) {
  
  model <- readRDS(x)  
  
  tibble::tibble(
    station = basename(x) %>% 
      stringr::str_sub(1, 8),
    depth = basename(x) %>% 
      stringr::str_sub(9, 10),
    params = summary(model[[1]])$coefficients %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="parameter") %>% 
      dplyr::select(parameter, value=Estimate) %>% 
      tidyr::pivot_wider(names_from=parameter) %>% 
      jsonlite::toJSON(digits=8),
    model = basename(x) %>% 
      stringr::str_split_1("_") %>% 
      magrittr::extract(2)
  ) %>% 
    dplyr::relocate(model, .before=params) 
}

raw_data_to_csv <- function(x) {
  model <- readRDS(x)  
  
  tibble::tibble(
    station = basename(x) %>% 
      stringr::str_sub(1, 8),
    depth = basename(x) %>% 
      stringr::str_sub(9, 10),
    data = model[[2]]
  ) %>% 
    tidyr::unnest(data) %>% 
    dplyr::rename(
      kpa=raw_kPa,
      vwc=raw_VWC
    )
}

process_all_to_table <- function(func, data_dir, out_dir, fname) {
  out <- list.files(data_dir, full.names = T, pattern = ".RDS") %>% 
    purrr::map(func)
  
  out %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(
      depth = dplyr::recode(
        depth, 
        "02"=-100,
        "04"=-100,
        "08"=-200,
        "20"=-500
      )
    ) %>%
    dplyr::distinct() %>% 
    readr::write_csv(file.path(out_dir, fname))
}

reformat_error_table <- function(f = "./data/Error_MSE/MSE", out_dir) {
  readr::read_csv(f, show_col_types = FALSE) %>% 
    dplyr::mutate(
      Depth = dplyr::recode(
        Depth, 
        "02"=-100,
        "04"=-100,
        "08"=-200,
        "20"=-500
      ),
      model = dplyr::case_when(
        model == "Fredlund-Xing Model" ~ "FX", 
        TRUE ~ "VG"
      )
    ) %>% 
    janitor::clean_names() %>% 
    dplyr::relocate(model, .before = mse) %>% 
    dplyr::distinct() %>% 
    readr::write_csv(file.path(out_dir, "model_error.csv"))
}

process_all_to_table(params_to_csv, "./data/RDS", "./data/for_db/", "parameters.csv")
process_all_to_table(raw_data_to_csv, "./data/RDS", "./data/for_db/", "vwc_kpa_observations.csv")
reformat_error_table("./data/Error_MSE/MSE", "./data/for_db")
