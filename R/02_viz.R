movi_data <- read_csv("results/movi_data.csv")

viz_data <- data.frame(country = rep(country_tags, each=length(years)),
                       year = rep(years, each=length(country_tags)),
                       movi = unlist(movi_data))
