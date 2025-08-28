#install.packages(magrittr)
#install.packages(ggplot2)
install.packages("plotly")

library(magrittr)
library(ggplot2)

# Need to be on campus or on the VPN to connect to DB.
con <-
  DBI::dbConnect(RPostgres::Postgres(),
                 host = "fcfc-mesonet-db2.cfc.umt.edu",
                 dbname = "mesonet",
                 user = "mesonet",
                 password = "Climate2019!"
  )

# This command will take a few seconds depending on how long of a record you want
dat <- dplyr::tbl(con, RPostgres::Id(schema = "api", table = "ln")) %>%
  dplyr::filter(
    station %in% c("acemedic"),
    # can change "soil_temp" to "soil_vwc", etc.
    stringr::str_detect(element, "soil_temp"),
    datetime >= "2025-05-01"
  ) %>%
  dplyr::collect()

dat %>% 
  ggplot(aes(x=datetime, y=value, color=id)) + 
  geom_line() + 
  facet_wrap(~element)
