# reading the data

# Raw species data
raw_bfly_Flona <- read.csv(here::here("data/Processed/data_bfly.csv"), sep = ";")
str(raw_bfly_Flona)
raw_bfly_Flona$Code <- paste(substr(raw_bfly_Flona$Strata, 1,1), raw_bfly_Flona$SU, sep = "")

# Raw environmental data
raw_env_Flona <- read.csv(here::here('data/Processed/environmental_data.csv'), sep = ";")
str(raw_env_Flona)
raw_env_Flona$Code <- paste(substr(raw_env_Flona$Strata, 1,1), raw_env_Flona$SU, sep = "")

# merging both data
library(dplyr)
tmp <- subset(raw_env_Flona, Period == "Day")
env_grouped <- data.frame(Tmean.day = round(tapply(tmp$Temperature, tmp$Code, mean, na.rm = T), 3),
                          Tsd.day = round(tapply(tmp$Temperature, tmp$Code, sd, na.rm = T), 3),
                          Hmean.day = round(tapply(tmp$Humidity, tmp$Code, mean, na.rm = T), 3),
                          Hsd.day = round(tapply(tmp$Humidity, tmp$Code, sd, na.rm = T), 3))

tmp <- subset(raw_env_Flona, Period == "Night")
env_grouped.n <- data.frame(Tmean.night = round(tapply(tmp$Temperature, tmp$Code, mean, na.rm = T), 3),
                            Tsd.night = round(tapply(tmp$Temperature, tmp$Code, sd, na.rm = T), 3),
                            Hmean.night = round(tapply(tmp$Humidity, tmp$Code, mean, na.rm = T), 3),
                            Hsd.night = round(tapply(tmp$Humidity, tmp$Code, sd, na.rm = T), 3))


env_grouped <- cbind(env_grouped, env_grouped.n)
env_grouped$Code <- rownames(env_grouped)

data_flona <- left_join(raw_bfly_Flona, env_grouped, by = "Code" )
data_flona %>% head

# extracting the values for richness and abundance
data_flona$Abundance <- rowSums(data_flona[,4:41])

# extracting the richness by day
data_flona$Richness <- rowSums(vegan::decostand(data_flona[,4:41], method = "pa"))

data_flona %>% head

data <- data_flona %>% na.omit()
data %>% head()
saveRDS(data, here::here("output/unscaled_data.rds"))

# scaling the predictor variables
data$Tmean.day <- scale(data$Tmean.day)
data$Tmean.night <- scale(data$Tmean.night)
data$Hmean.day <- scale(data$Hmean.day)
data$Hmean.night <- scale(data$Hmean.night)
data$Tsd.day <- scale(data$Tsd.day)
data$Tsd.night <- scale(data$Tsd.night)
data$Hsd.day <- scale(data$Hsd.day)
data$Hsd.night <- scale(data$Hsd.night)

# evaluating the colinearity graphically
source(here::here("R/functions/panel_functions.R"))

pairs(data[,43:50], upper.panel=panel.cor, diag.panel=panel.hist)

# saving the data
colnames(data)
data_models <- data[,c("Strata", "SU", "Tmean.day", "Tsd.day", 
                       "Hmean.day", "Hsd.day", "Tmean.night", "Tsd.night",
                       "Hmean.night", "Hsd.night", "Abundance", "Richness")]

data_comm <- data[,c(1:50)]
data_comm %>% head()

saveRDS(object = data_models, file = here::here("output/data_models_std.rds"))
saveRDS(data_comm, here::here("output/community_env_data.rds"))
