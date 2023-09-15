movi_data <- read_csv("results/movi_data.csv")

viz_data <- data.frame(country = rep(country_tags, each=length(years)),
                       year = rep(years, length(country_tags)),
                       movi = unlist(movi_data))

ggplot(data = viz_data, aes(x=year, y = movi, group = country, colour=country)) + geom_line(size = 1) 
#+ facet_wrap(~country,scales = "free_y", ncol = 2)

ggplot(viz_data, aes(fill=country, y=movi, x=year)) + 
  geom_bar(position="stack", stat="identity")


movi_averages <- data.frame(country = country_tags,
                          average_global = rep(NA,length(country_tags)),
                          average_10 = rep(NA,length(country_tags)),
                          average_5 = rep(NA,length(country_tags)))
i <- 1
for (country_i in country_tags){
  temp_data <- viz_data %>% filter(country == country_i)
  temp_avr_global <- mean(temp_data$movi)
  temp_data <- viz_data %>% filter(country == country_i & year > 2012)
  temp_avr_10 <- mean(temp_data$movi)
  temp_data <- viz_data %>% filter(country == country_i & year > 2017)
  temp_avr_5 <- mean(temp_data$movi)
  
  movi_averages[names(country_tags)[i],"average_global"] <- temp_avr_global
  movi_averages[names(country_tags)[i],"average_10"] <- temp_avr_10
  movi_averages[names(country_tags)[i],"average_5"] <- temp_avr_5
  i <- i+1
}

viz_averages <- viz_data %>% mutate(
  movi_ave_global = rep(movi_averages$average_global,each = length(years)),
  movi_ave_10 = rep(movi_averages$average_10,each = length(years)),
  movi_ave_5 = rep(movi_averages$average_5,each = length(years)),
  change_global = movi/movi_ave_global - 1,
  change_10 = movi/movi_ave_10 - 1,
  change_5 = movi/movi_ave_5 - 1
)

ggplot(viz_averages, aes(fill=country, y=change_5, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  geom_hline(yintercept = 0) +  facet_wrap(~country,scales = "free_y", ncol = 2)

write.csv(viz_averages, file = 'results/movi_data_averages.csv', row.names = F)
