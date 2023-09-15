movi_data <- read_csv("results/movi_data.csv")

viz_data <- data.frame(country = rep(country_tags, each=length(years)),
                       year = rep(years, length(country_tags)),
                       movi = unlist(movi_data))

ggplot(data = viz_data, aes(x=year, y = movi, group = country, colour=country)) + geom_line(size = 2)
