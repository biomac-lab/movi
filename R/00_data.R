library(readr)
library(tidyverse)

country_tags <- c("ARG" = "Argentina", "BRA" = "Brazil", "BOL" = "Bolivia",
                  "CHL" = "Chile", "COL" = "Colombia", "CRI" = "Costa Rica",
                  "ECU" ="Ecuador", "GTM" = "Guatemala", "HND" = "Honduras",
                  "MEX" = "Mexico", "NIC" = "Nicaragua", "PAN" = "Panama",
                  "PER" = "Peru", "PRY" = "Paraguay", "SLV" = "El Salvador",
                  "URY" = "Uruguay", "VEN" = "Venezuela")

years <- seq(2000,2022)


## Climate data

t_data_path <- c("data/temperature_data/tas_timeseries_annual_era_1970-2020_")
h_data_path <- c("data/humidity_data/hurs_timeseries_annual_era_1970-2020_")
p_data_path <- c("data/precipitation_data/pr_timeseries_annual_era_1970-2020_")

climate_data <- data.frame("year" = rep(years, length(country_tags)),
                           "country" = rep(country_tags, each = length(years)), 
                           "iso3" = rep(names(country_tags), each = length(years)),
                           "temperature" = rep(NA, length(country_tags)),
                           "humidity" = rep(NA, length(country_tags)),
                           "precipitation" = rep(NA, length(country_tags)))

# Temperature data
t_data <- c()
for (tag in names(country_tags)){
  file <- paste0(t_data_path,tag,".csv")
  orig_data <- as.data.frame(read_csv(file, na = "NA", skip = 1))
  orig_data <- orig_data %>% mutate(year = orig_data[[1]])
  orig_data <- orig_data %>% filter( year %in% years) %>% select(2,"year")
  reps <- rep(last(orig_data[[1]]),last(years)-last(orig_data[[2]]))
  imported_data <- c(orig_data[[1]],reps)
  t_data <- c(t_data,imported_data)
}

# Humidity data
h_data <- c()
for (tag in names(country_tags)){
  file <- paste0(h_data_path,tag,".csv")
  orig_data <- as.data.frame(read_csv(file, na = "NA", skip = 1))
  orig_data <- orig_data %>% mutate(year = orig_data[[1]])
  orig_data <- orig_data %>% filter( year %in% years) %>% select(2,"year")
  reps <- rep(last(orig_data[[1]]),last(years)-last(orig_data[[2]]))
  imported_data <- c(orig_data[[1]],reps)
  h_data <- c(h_data,imported_data)
}

# Precipitation data
p_data <- c()
for (tag in names(country_tags)){
  file <- paste0(p_data_path,tag,".csv")
  orig_data <- as.data.frame(read_csv(file, na = "NA", skip = 1))
  orig_data <- orig_data %>% mutate(year = orig_data[[1]])
  orig_data <- orig_data %>% filter( year %in% years) %>% select(2,"year")
  reps <- rep(last(orig_data[[1]]),last(years)-last(orig_data[[2]]))
  imported_data <- c(orig_data[[1]],reps)
  p_data <- c(p_data,imported_data)
}

climate_data$temperature <- t_data
climate_data$humidity <- h_data
climate_data$precipitation <- p_data


## HAQ Data
haq_orig_data <- as.data.frame(read_delim("data/haq_data.csv",
                                          delim = ";", escape_double = FALSE))
haq_orig_data <- haq_orig_data %>% filter(year %in% years & iso3 %in% names(country_tags))

haq_data <- data.frame()

for (tag in names(country_tags)){
  orig_data <- haq_orig_data %>% filter(iso3 == tag)
  reps <- as.data.frame(rep(last(orig_data),last(years)-max(orig_data$year)))
  reps$year <- seq(max(orig_data$year)+1,last(years))
  haq_data <- rbind(haq_data,orig_data,reps)
}

## WASH Data

wash_orig_data <- as.data.frame(read_delim("data/wash_data.csv", 
                                           delim = ";", escape_double = FALSE, na = "NA"))

service_lvl_tags <- c("Surface water" = 1, "Unimproved" = 2,
                      "Limited service" = 3, "At least basic" = 4,
                      "Basic service" = 5, "Safely managed service" = 6)

coverages <- c()
for (c in names(country_tags)){
  for (y in years){
    # sum_cov <- sum(wash_orig_data$coverage[which(wash_orig_data$service_level %in% 
    #                                                names(service_lvl_tags[1:3]) &
    #                                                wash_orig_data$iso3 == c &
    #                                                wash_orig_data$year== y)])
    # TODO: corregir suma
    print(y)
    print(sum_cov)
    coverages <- c(coverages, sum_cov)
  }
}

wash_data <- data.frame(iso3 = rep(unique(wash_orig_data$iso3), each = length(years)),
                        country = rep(unique(wash_orig_data$country), each = length(years)),
                        year = rep(years, length(country_tags)),
                        coverage = coverages)

## Altitudinal demographic data

alt_orig_data <- as.data.frame(read_delim("data/altitudinal_demog_data.csv", 
                            delim = ";", escape_double = FALSE, na = "NA"))
alt_orig_data$`<500m` <- as.numeric(alt_orig_data$`<500m`)

alt_data <- alt_orig_data %>% filter(country %in% country_tags) %>% 
  mutate(
    n_below2000m = `<500m` + `500-999m` + `1000-1499m` + `1500-1999m`,
    p_below2000m = `%<500m` + `%500-999m` + `%1000-1499m` + `%1500-1999m`
  )


## Urban data

urban_orig_data <- read_delim("data/urban_population_data.csv", 
                                    delim = ";", escape_double = FALSE, na = "NA")

urban_orig_data <- urban_orig_data %>% filter(year %in% years & iso3 %in% names(country_tags))

urban_data <- data.frame()

for (tag in names(country_tags)){
  orig_data <- urban_orig_data %>% filter(iso3 == tag)
  reps <- as.data.frame(rep(last(orig_data),last(years)-max(orig_data$year)))
  reps$year <- seq(max(orig_data$year)+1,last(years))
  urban_data <- rbind(urban_data,orig_data,reps)
}

#rm(list=(ls()[!(ls() %in% c("climate_data","haq_data","urban_data","wash_data","alt_data",
                            "country_tags", "years", "service_lvl_tags"))]))
