library(readr)
library(tidyverse)
library(readxl)

country_tags <- c("ARG" = "Argentina", "BOL" = "Bolivia", "BRA" = "Brazil", 
                  "CHL" = "Chile", "COL" = "Colombia", "CRI" = "Costa Rica",
                  "ECU" ="Ecuador", "SLV" = "El Salvador", "GTM" = "Guatemala",
                  "HND" = "Honduras", "MEX" = "Mexico", "NIC" = "Nicaragua",
                  "PAN" = "Panama", "PRY" = "Paraguay", "PER" = "Peru",  
                  "URY" = "Uruguay", "VEN" = "Venezuela")

years <- seq(2000,2023)


## Climate data

climate_data <- data.frame("year" = rep(years, length(country_tags)),
                           "country" = rep(country_tags, each = length(years)),
                           "iso3" = rep(names(country_tags), each = length(years)),
                           "temperature" = rep(NA, length(country_tags)),
                           "humidity" = rep(NA, length(country_tags)),
                           "precipitation" = rep(NA, length(country_tags)))

# Temperature data

tas_1950_2023 <- read_excel(paste0(getwd(),"/data/temperature_data/era5_tas_1950-2023.xlsx"))
colnames(tas_1950_2023) <- c("code", "name", as.character(seq(1950,2023)))

t_data <- c()
for (tag in names(country_tags)){
  tas_country <- subset(tas_1950_2023, code == tag)
  orig_data <- as.vector(select(tas_country, as.character(years)))
  reps <- rep(last(orig_data),last(years)-as.numeric(last(names(orig_data))))
  complete_data <- c(orig_data,reps)
  t_data <- c(t_data,complete_data)
}

# Humidity data
hurs_1950_2023 <- read_excel(paste0(getwd(),"/data/humidity_data/era5_hurs_1950-2023.xlsx"))
colnames(hurs_1950_2023) <- c("code", "name", as.character(seq(1950,2023)))

h_data <- c()
for (tag in names(country_tags)){
  h_country <- subset(hurs_1950_2023, code == tag)
  orig_data <- as.vector(select(h_country, as.character(years)))
  reps <- rep(last(orig_data),last(years)-as.numeric(last(names(orig_data))))
  complete_data <- c(orig_data,reps)
  h_data <- c(h_data,complete_data)
}

# Precipitation data
pr_1950_2023 <- read_excel(paste0(getwd(),"/data/precipitation_data/era5_pr_1950-2023.xlsx"))
colnames(pr_1950_2023) <- c("code", "name", as.character(seq(1950,2023)))

p_data <- c()
for (tag in names(country_tags)){
  p_country <- subset(pr_1950_2023, code == tag)
  orig_data <- as.vector(select(p_country, as.character(years)))
  reps <- rep(last(orig_data),last(years)-as.numeric(last(names(orig_data))))
  complete_data <- c(orig_data,reps)
  p_data <- c(p_data,complete_data)
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
  orig_data <- as.data.frame(haq_orig_data %>% filter(iso3 == tag))
  reps <- do.call(rbind, c(replicate(last(years)-max(orig_data$year), orig_data[nrow(orig_data),], simplify = FALSE)))
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
    country_cov <- wash_orig_data %>% filter(iso3 == c, wash_orig_data$year== y,  
                                             wash_orig_data$service_level %in% names(service_lvl_tags[1:3])) %>%
                                              select(coverage)
    while(nrow(country_cov)==0){
      y <- y-1
      country_cov <- wash_orig_data %>% filter(iso3 == c, wash_orig_data$year== y,  
                                               wash_orig_data$service_level %in% names(service_lvl_tags[1:3])) %>%
        select(coverage)
    }
    coverages <- c(coverages, sum(country_cov))
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
  orig_data <- as.data.frame(urban_orig_data %>% filter(iso3 == tag))
  reps <- do.call(rbind, c(replicate(last(years)-max(orig_data$year), orig_data[nrow(orig_data),], simplify = FALSE)))
  reps$year <- seq(max(orig_data$year)+1,last(years))
  urban_data <- rbind(urban_data,orig_data,reps)
}
