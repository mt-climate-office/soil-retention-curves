library(magrittr)
source("./R/fit_curve.R")

recode_depth <- function(dat) {
  dat %>% 
    dplyr::mutate(
      depth = stringr::str_replace(depth, "cm", "") %>% 
        as.numeric() %>% 
        magrittr::multiply_by(-10)
    )
}

params_to_csv <- function(model, station, depth) {
  
  tibble::tibble(
    station = station,
    depth = depth,
    params = summary(model$model)$coefficients %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var="parameter") %>% 
      dplyr::mutate(
        parameter = stringr::str_replace(parameter, "O_", "")
      ) %>%
      dplyr::select(parameter, value=Estimate) %>% 
      tidyr::pivot_wider(names_from=parameter) %>% 
      jsonlite::toJSON(digits=8) %>% 
      list(),
    model = "FX"
  ) %>% 
    dplyr::relocate(model, .before=params) %>%
    recode_depth()
}

raw_data_to_csv <- function(model, station, depth) {

  tibble::tibble(
    station = station,
    depth = depth,
    data = model$raw_data
  ) %>% 
    tidyr::unnest(data) %>% 
    dplyr::rename(
      kpa=raw_kPa,
      vwc=raw_VWC
    ) %>% 
    recode_depth()
}

# process_all_to_table <- function(func, data_dir, out_dir, fname) {
#   out <- list.files(data_dir, full.names = T, pattern = ".RDS") %>% 
#     purrr::map(func)
#   
#   out %>% 
#     dplyr::bind_rows() %>%
#     dplyr::mutate(
#       depth = dplyr::recode(
#         depth, 
#         "02"=-100,
#         "04"=-100,
#         "08"=-200,
#         "20"=-500
#       )
#     ) %>%
#     dplyr::distinct() %>% 
#     readr::write_csv(file.path(out_dir, fname))
# }

reformat_error_table <- function(model, station, depth) {
  model$error %>% 
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "model", 
      values_to = "mse"
    ) %>% 
    dplyr::filter(model != "Kosugi_RSE") %>%
    dplyr::mutate(
      station = station, 
      depth = depth,
      model = dplyr::case_when(
        model == "van_RSE" ~ "VG",
        TRUE ~ "FX"
      ),
      mse = mse * mse
    ) %>%
    recode_depth()
}


process_zipped_data <- function(data_dir) {
  
  models <- list.files(data_dir, full.names = T, pattern = ".zip") %>%
    purrr::map(fit_all_depths) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rowwise()
  
  params <- models %>% 
    dplyr::mutate(tmp = list(params_to_csv(model, station, depth))) %>% 
    dplyr::select(tmp) %>% 
    tidyr::unnest(cols = tmp) %>% 
    tidyr::unnest(cols = params)
  
  raw <- models %>% 
    dplyr::mutate(tmp = list(raw_data_to_csv(model, station, depth))) %>% 
    dplyr::select(tmp) %>% 
    tidyr::unnest(cols = tmp)
  
  error <- models %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(tmp = list(reformat_error_table(model, station, depth))) %>% 
    dplyr::select(tmp) %>% 
    tidyr::unnest(cols = tmp)
  
  porosity <- models %>%
    dplyr::select(station, depth, porosity) %>%
    recode_depth()
  
  list(
    "params" = params,
    "raw" = raw,
    "error" = error,
    "porosity" = porosity
  )
}

write_to_db <- function(tables) {
  con <-
    DBI::dbConnect(RPostgres::Postgres(),
                   host = Sys.getenv("MESONET_HOSTNAME"),
                   dbname = Sys.getenv("MESONET_DBNAME"),
                   user = Sys.getenv("MESONET_USER"),
                   password = Sys.getenv("MESONET_PASSWORD")
    )

  tryCatch({
    DBI::dbAppendTable(conn = con,
                       name = DBI::Id(schema = "soil",
                                      table = "model_error"),
                       value = tables$error)
  }, error = function(e) {
    if (grepl("unique constraint", e$message, ignore.case = TRUE)) {
      print("Records already exist in error table. Deleting original and uploading updated records.")
      station_names <- unique(tables$error$station)
      stations_str <- paste(shQuote(station_names), collapse = ", ")
      query <- sprintf("DELETE FROM soil.model_error WHERE station IN (%s)", stations_str)
      DBI::dbExecute(con, query)
      
      DBI::dbAppendTable(conn = con,
                         name = DBI::Id(schema = "soil",
                                        table = "model_error"),
                         value = tables$error)
    } else {
      stop(e)
    } 
  })

  
  tryCatch({
    DBI::dbAppendTable(conn = con,
                       name = DBI::Id(schema = "soil",
                                      table = "parameters"),
                       value = tables$params)
  }, error = function(e) {
    if (grepl("unique constraint", e$message, ignore.case = TRUE)) {
      print("Records already exist in parameters table. Deleting original and uploading updated records.")
      station_names <- unique(tables$params$station)
      stations_str <- paste(shQuote(station_names), collapse = ", ")
      query <- sprintf("DELETE FROM soil.parameters WHERE station IN (%s)", stations_str)
      DBI::dbExecute(con, query)
      
      DBI::dbAppendTable(conn = con,
                         name = DBI::Id(schema = "soil",
                                        table = "parameters"),
                         value = tables$params)
    } else {
      stop(e)
    } 
  })
  
  tryCatch({
    DBI::dbAppendTable(conn = con,
                       name = DBI::Id(schema = "soil",
                                      table = "porosity"),
                       value = tables$porosity)
  }, error = function(e) {
    if (grepl("unique constraint", e$message, ignore.case = TRUE)) {
      print("Records already exist in porosity table. Deleting original and uploading updated records.")
      station_names <- unique(tables$porosity$station)
      stations_str <- paste(shQuote(station_names), collapse = ", ")
      query <- sprintf("DELETE FROM soil.porosity WHERE station IN (%s)", stations_str)
      DBI::dbExecute(con, query)
      
      DBI::dbAppendTable(conn = con,
                         name = DBI::Id(schema = "soil",
                                        table = "porosity"),
                         value = tables$porosity)
    } else {
      stop(e)
    } 
  })


  station_names <- unique(tables$raw$station)
  stations_str <- paste(shQuote(station_names), collapse = ", ")
  query <- sprintf("DELETE FROM soil.raw_soil_data WHERE station IN (%s)", stations_str)
  DBI::dbExecute(con, query)
  
  DBI::dbAppendTable(conn = con,
                     name = DBI::Id(schema = "soil",
                                    table = "raw_soil_data"),
                     value = tables$raw)

  
  
  DBI::dbDisconnect(con)
}


# tables <- process_zipped_data("./data/completed_uploads")
the_whole_thing <- function(data_dir) {
  process_zipped_data(data_dir) %>%
    write_to_db()
  
  tibble::tibble(
    from = list.files(data_dir, full.names = T)
  ) %>% 
    dplyr::mutate(
      to = stringr::str_replace(from,"zipped_data", "completed_uploads")
    ) %>%
    purrr::pmap(file.rename)
}

the_whole_thing("./data/zipped_data")

upload_old_porosity <- function() {
  out <- list.files("./data/Raw_Data/", full.names = T) %>%
    purrr::map(function(x) {
      print(x)
      
      porosity = get_porosity(x)
      
      # depth = str_extract(x, "(?<=_)(02|04|08|20|36)(?=_)") %>%
      depth = x %>%
        basename() %>% 
        tools::file_path_sans_ext() %>% 
        stringr::str_sub(9) %>%
        dplyr::case_match(
          "02" ~ "10cm",
          "04" ~ "10cm",
          "08" ~ "20cm",
          "20" ~ "50cm", 
          "36" ~ "91cm"
        )
      
      
      station <- x %>%
        basename() %>% 
        tools::file_path_sans_ext() %>%
        stringr::str_sub(end=-3) 
      
      if (is.na(depth)) {
        stop(glue::glue("No depth for {station}"))
      }
      
      
      tibble::tibble(
        station = station,
        depth =depth,
        porosity=porosity
      )
    }) %>%
    dplyr::bind_rows() %>%
    recode_depth()
  
  
  tryCatch({
    DBI::dbAppendTable(conn = con,
                       name = DBI::Id(schema = "soil",
                                      table = "porosity"),
                       value = out)
  }, error = function(e) {
    if (grepl("unique constraint", e$message, ignore.case = TRUE)) {
      print("Records already exist in porosity table. Deleting original and uploading updated records.")
      station_names <- unique(out$station)
      stations_str <- paste(shQuote(station_names), collapse = ", ")
      query <- sprintf("DELETE FROM soil.porosity WHERE station IN (%s)", stations_str)
      DBI::dbExecute(con, query)
      
      DBI::dbAppendTable(conn = con,
                         name = DBI::Id(schema = "soil",
                                        table = "porosity"),
                         value = out)
    } else {
      stop(e)
    } 
  })
  
  
}
