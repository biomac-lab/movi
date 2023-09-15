source(file = "R/00_data.R")
source(file = "R/MVSE_package.R")

require("pbapply")
require("scales")
require("genlasso")
require("readr")


## Data exploration

ggplot(data = urban_data, aes(x=year, y = urban_population, group = country, colour=country)) + geom_line()
ggplot(data = haq_data, aes(x=year, y = hca, group = country, colour=country)) + geom_line()
ggplot(data = wash_data, aes(x=year, y = coverage, group = country, colour=country)) + geom_line()


## Index P estimation

indexP_data <- data.frame(matrix(ncol = length(country_tags),
                               nrow = length(years)))
colnames(indexP_data) <- country_tags[]

for (country_i in country_tags){
  
  print(country_i)
  
  ClimateSeries <- climate_data %>% filter(country == country_i)
  ClimateSeries <- data.frame(T = ClimateSeries$temperature,
                              H = ClimateSeries$humidity,
                              year = ClimateSeries$year,
                              date = as.Date.character(paste0(as.character(ClimateSeries$year),"-01","-01"), format = "%Y-%m-%d"),
                              R = ClimateSeries$precipitation)
  
  climate_file_path <- paste0("data/climate_ready_",country_i,".csv")
  write.csv(ClimateSeries, file = climate_file_path, row.names = F)
  
  setEmpiricalClimateSeries(climate_file_path)
  setOutputFilePathAndTag(paste0('results/movi_wP_',country_i))
  
  plotClimate()
  
  setMosqLifeExpPrior(pmean=12, psd=2, pdist='gamma')  
  setMosqIncPerPrior(pmean=7, psd=2, pdist='gamma')  
  setMosqBitingPrior(pmean=0.25, psd=0.01, pdist='gamma')  
  setHumanLifeExpPrior(pmean=71.1, psd=2, pdist='gamma')
  setHumanIncPerPrior(pmean=5.8, psd=1, pdist='gamma')
  setHumanInfPerPrior(pmean=5.9, psd=1, pdist='gamma')
  setHumanMosqTransProbPrior(pmean=0.5, psd=0.01, pdist='gamma')
  
  estimateEcoCoefficients(nMCMC=100000, bMCMC=0.5, cRho=1, cEta=1, gauJump=0.75)
  
  simulateEmpiricalIndexP(nSample=1000, smoothing=c(7,15,30,60))
  
  exportEmpiricalIndexP()
  plotEmpiricalIndexP(outfilename=paste0('debug_indexP',country_i))
  
  indexP_country <- as.data.frame(read_csv(paste0("results/movi_wP_",country_i,
                                                  ".estimated_indexP.csv")))
  indexP_country <- as.numeric(indexP_country$indexP)
  
  indexP_data <- indexP_data %>% mutate(
    !!country_i := indexP_country
  )
  
}

indexP_data_sca <- indexP_data/sapply(indexP_data, max, na.rm = TRUE)
indexP_data_sca <- indexP_data_sca %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))

movi_data <- data.frame(matrix(ncol = length(country_tags),
                               nrow = length(years)))
colnames(movi_data) <- country_tags[]

for (country_i in country_tags){
  
  indexP_country <- as.numeric(indexP_data_sca[,country_i])
  
  urban_country <- urban_data %>% filter(country == country_i)
  urban_country <- as.numeric(urban_country$urban_population)/100
  
  wash_country <- wash_data %>% filter(country == country_i)
  wash_country <- as.numeric(wash_country$coverage)
  
  haq_country <- haq_data %>% filter(country == country_i) 
  haq_country <- as.numeric(haq_country$hca)/100
  
  alt_country <- alt_data %>% filter(country == country_i)
  alt_country <- rep(as.numeric(alt_country$p_below2000m), length(years))/100
  
  movi_country <- 100*indexP_country*(urban_country*wash_country*alt_country)/haq_country
  
  movi_data <- movi_data %>% mutate(
    !!country_i := movi_country
  )
  
}
write.csv(movi_data, file = 'results/movi_data.csv', row.names = F)

