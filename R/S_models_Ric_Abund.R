

# Required packages -------------------------------------------------------
library(ggeffects)
library(ggplot2)
library(viridis)

# reading the data
abn.nb <- readRDS(here::here("output/output_model_abn.rds"))
data_models <- readRDS(here::here("output/data_models_std.rds"))

summary(abn.nb) # model for abundance

# predicting the values
pred.abn.tmean <- ggpredict(abn.nb, c("Tmean.day", "Strata"), type = "re")
pred.abn.tmean$group

plot(pred.abn.tmean)

p.tmean <- ggplot(pred.abn.tmean, aes(x = x, y = predicted, colour = group)) +
  geom_smooth(method = "loess") + 
  geom_point(data = data_models, aes(x = Tmean.day, y = Abundance, colour = Strata), 
                                   alpha = .5) +
 # geom_ribbon(aes(ymin= conf.low, ymax= conf.high),
 #                      alpha = .1, linetype = 0) + 
  scale_y_sqrt() +   scale_color_manual(values = c("firebrick3", "palegreen3" )) +
  labs(x = "Mean Temperature", y = "Predicted Abundance", tag = "a)")
p.tmean
# quando a temperature média do aumenta, a tendencia é aumantar a abundancia das espécies,
# independente do estrato amostrado

pred.abn.tsdn <- ggpredict(abn.nb, c("Tsd.night", "Strata"), type = "re")
plot(pred.abn.tsdn)

p.tsd <- ggplot(pred.abn.tsdn, aes(x = x, y = predicted, colour = group)) +
  geom_smooth() + 
  geom_point(data = data_models, aes(x = Tsd.night, y = Abundance, colour = Strata), 
             alpha = .5) +
 #geom_ribbon(aes(ymin= conf.low, ymax= conf.high,  fill = group),
         #           alpha = .1, linetype = 0) + 
  scale_y_sqrt() + 
  scale_color_manual(values = c("firebrick3", "palegreen3" )) + #scale_fill_viridis_d(option = "B")  +
  labs(x = "Temperature variation", y = "Predicted Abundance", colour = "Strata", 
       fill = "Strata", tag = "b)")
p.tsd

# a variação na temperatura da noite afeta a abudancia das espécies dependendo do estrtao onde elas
# foram amostradas, nesse caso, menores variações de temperatura na noite levam a maiores abundancias 
# no subosque, enquanto altas variações de temperatura levam a um aumento da abundancia das espécies no dossel

# interaction between Strata and Hmean.day
pred.abn.hmean <- ggpredict(abn.nb, c("Hmean.day", "Strata"), type = "re")
plot(pred.abn.hmean)

p.hmean <- ggplot(pred.abn.hmean, aes(x = x, y = predicted, colour = group)) +
  geom_smooth() + 
  geom_point(data = data_models, aes(x = Hmean.day, y = Abundance, colour = Strata),
             alpha = .5) +
  scale_y_sqrt() +
  scale_color_manual(values = c("firebrick3", "palegreen3" )) +
  labs(x = "Mean Humidity", y = "Predicted Abundance", colour = "Strata",
       fill = "Strata", tag = "c)")
p.hmean 

# a media de umidade afeta a abundancia das espécies dependendo do estrato em que elas foram 
# amostradas, baixas medias levam a menores abundancias no dossel, porém altas medias aumentam
# a abudancia no dossel. baixa umidade restringe a quantidade de individuos no dossel.

abn.models <- cowplot::plot_grid(p.tmean + theme(legend.position = "none"), p.tsd + theme(legend.position = "none"),
                   p.hmean + theme(legend.position = "none"))

cowplot::save_plot(here::here("output/Images/G_abundance_models.png"),
                   abn.models, base_height = 8, base_width = 10)

##############################################################################
##############################################################################
##############################################################################
# trying to extract confidence intervals for Hmean.day:Strata -------------

colnames(pred.abn.hmean)
colnames(pred.abn.tmean)
pred.abn.tmean$std.error

library(lme4)
library(ggplot2) # Plotting
data("Orthodont",package="nlme")

newdat <- expand.grid(
  Hmean.day = c(-1.2, -.6, 0, .6, 1.2),
  Strata =c("Canopy","Underestory"),
  Tmean.day = c(-1.2, -.6, 0, .6, 1.2),
  Tsd.night = c(-1.2, -.6, 0, .6, 1.2),
  Abundance = 0)

mm <- model.matrix(terms(abn.nb), newdat)
newdat$Abundance <- mm %*% fixef(abn.nb)

pvar1 <- diag(mm %*% tcrossprod(vcov(abn.nb),mm))
tvar1 <- pvar1+VarCorr(abn.nb)$SU[1]  ## must be adapted for more complex models

newdat <- data.frame(
  newdat
  , plo = newdat$Abundance-2*sqrt(pvar1) #fixed effects only
  , phi = newdat$Abundance+2*sqrt(pvar1)
  , tlo = newdat$Abundance-2*sqrt(tvar1) # fixed + random effects
  , thi = newdat$Abundance+2*sqrt(tvar1)
)

head(newdat)

#plot confidence
g0 <- ggplot(newdat, aes(x=Tmean.day, y=Abundance))+
  geom_line()
g0 + geom_ribbon(aes(ymin= plo, ymax= phi,  fill = Strata),
                 alpha = .1, linetype = 0) + scale_color_viridis_d() +
  scale_fill_viridis_d()
#plot prediction
g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
  labs(title="CI based on FE uncertainty + RE variance")

# manually calculate the confidence intervals
# +/- 1.96 * SE
obj <- summary(abn.nb)[["coefficients"]]

tmp <- subset(pred.abn.tsdn, group == "Canopy")
ci.up.c <- tmp$predicted + (2 * exp(obj[3,2])) 
ci.lo.c <- tmp$predicted - (2 * exp(obj[3,2]))

tmp <- subset(pred.abn.tsdn, group == "Understory")
ci.up.u  <- tmp$predicted + exp(1.96 * (obj[5,2]))
ci.lo.u  <- tmp$predicted - (1.96 * exp(obj[5,2]))

t <- cbind(conf.low = c(ci.lo.c, ci.lo.u), vec = gl(14, 1))
conf.low <- t[order(t[,2]),1]

pred.abn.tsdn
pred.abn.tsdn$conf.low.1 <- conf.low

t <- cbind(conf.low = c(ci.up.c, ci.up.u), vec = gl(14, 1))
conf.high <- t[order(t[,2]),1]

pred.abn.tsdn$conf.high.1 <- conf.high


p.test <- ggplot(pred.abn.tsdn, aes(x = x, y = predicted, colour = group)) +
  geom_line() + 
  geom_point(data = data_models, aes(x = Tsd.night, y = Abundance, colour = Strata)) +
  geom_ribbon(aes(ymin= conf.low.1, ymax= conf.high.1,  fill = group),
              alpha = .1, linetype = 0) + scale_color_viridis_d() +
  scale_fill_viridis_d()
p.test




  geom_ribbon(aes(ymin= conf.low, ymax= conf.high,  fill = cat.hum),
              alpha = .1, linetype = 0)  + guides(colour = F, fill = F) +
  #scale_color_viridis_d() +
  #scale_fill_viridis_d()+
  labs(x= "Temperature", y = "Predicted Richness", 
       colour = "Humidity", fill = "Humidity", tag = "a)")
plot.ric
