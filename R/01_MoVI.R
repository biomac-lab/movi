source(file = "R/00_data.R")
source(file = "R/MVSE_package.R")

require("pbapply")
require("scales")
require("genlasso")

country_i <- country_tags[5]

ClimateSeries <- climate_data %>% filter(country == country_i)
ClimateSeries <- data.frame(T = ClimateSeries$temperature,
                            H = ClimateSeries$humidity,
                            year = ClimateSeries$year,
                            date = as.Date.character(paste0(as.character(ClimateSeries$year),"-01","-01"), format = "%Y-%m-%d"),
                            R = ClimateSeries$precipitation)

write.csv(ClimateSeries, file = "data/climate_ready.csv", row.names = F)

setEmpiricalClimateSeries('data/climate_ready.csv')
setOutputFilePathAndTag('results/movi_wP')

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
plotEmpiricalIndexP(outfilename='debug_indexP')

COL_indexP <- as.data.frame(read_csv("results/movi_wP.estimated_indexP.csv"))
COL_indexP <- COL_indexP$indexP

COL_wash <- wash_data %>% filter(iso3 == "COL") %>% select(coverage)

COL_HAQ <- haq_data %>% filter(iso3 == "COL") %>% select(hca)

movi_col <- 1000*COL_indexP*COL_wash/COL_HAQ
