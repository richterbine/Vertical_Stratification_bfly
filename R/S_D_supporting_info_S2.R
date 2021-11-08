# script to evaluate the diferences in temeprature and humidity throughout day

data_models <- readRDS(here::here("output/data_models_std.rds"))
colnames(data_models)

df.data <- data.frame(Tvalue = c(data_models$Tmean.day, data_models$Tsd.day, 
                                data_models$Tmean.night, data_models$Tsd.night),
                      Hvalue = c(data_models$Hmean.day, data_models$Hsd.day, 
                                 data_models$Hmean.night, data_models$Hsd.night),
                      type = c(rep("Mean", nrow(data_models)), rep("SD", nrow(data_models))),
                      period = c(rep("Day", 2*nrow(data_models)), rep("Night", 2*nrow(data_models))),
                      Strata = rep(data_models$Strata, 2))

p.tvalue <- ggplot(data = df.data, aes(x = period, y = Tvalue, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = .5) + labs(x = NULL, tag = "a)", y = "Temperature") +
  scale_colour_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  facet_wrap(~type)
p.tvalue

p.hvalue <- ggplot(data = df.data, aes(x = period, y = Hvalue, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = .5) + labs(x = NULL, tag = "b)", y = "Humidity") +
  scale_colour_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
    facet_wrap(~type)

p.periods <- cowplot::plot_grid(p.tvalue + theme(legend.position = "none"), 
                   p.hvalue + theme(legend.position = "none"),
                   nrow = 2)

cowplot::save_plot(filename = here::here("output/Images/S2_temp_humi.png"), plot = p.periods,
                   base_height = 6, base_width = 8)

# evaluating the difference betwwen strata for climatic conditions
summary.tb <- list()

for (i in 3:10) {
  lm.tmp <- lm(data_models[,i] ~ Strata, data_models)
  summary.tb[[i]] <- summary(lm.tmp)$coefficients
}

df.summary <- as.data.frame(rbind(summary.tb[[3]], summary.tb[[4]],
                                  summary.tb[[5]], summary.tb[[6]],
                                  summary.tb[[7]], summary.tb[[8]],
                                  summary.tb[[9]], summary.tb[[10]]))
df.summary$type <- rep(colnames(data_models[,3:10]), each = 2)
write.csv(x = df.summary, file = here::here("output/env_var.csv"))
