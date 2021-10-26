
# Required packages -------------------------------------------------------
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(tidyverse)

# Evaluating the effect of environmental factors in Richness and Abundance --------
# read files
data_models <- readRDS(here::here("output/data_models_std.rds"))
data_models %>% str
sum(is.na(data_models))

###########################################################################

# For Richness data -------------------------------------------------------

# Null models
ric.pois.null <- glmer(Richness ~ 1 + (1|SU), family = poisson, data = data_models)

# create a formula allowing triple interactions
colnames(data_models)

ric.pois <- glmer(Richness ~ Strata + Tmean.day + Hmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(ric.pois)

# removing the Hmean.day
ric.pois <- glmer(Richness ~ Strata + Tmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(ric.pois)

# removing the Strata
ric.pois <- glmer(Richness ~ Tmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(ric.pois) # all <3

summary(ric.pois)
deviance(ric.pois)/df.residual(ric.pois)

r.squaredGLMM(ric.pois)
model.sel(ric.pois, ric.pois.null)

################ ABUNDANCE ##############################
# Null models
abn.pois.null <- glmer(Abundance ~ 1 + (1|SU), family = poisson, data = data_models)

abn.nb.null <- glmer.nb(Abundance ~ 1 + (1|SU), data = data_models)

# create a formula allowing triple interactions
colnames(data_models)
abn.pois <- glmer(Abundance ~ Strata + Tmean.day + Hmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(abn.pois)

# removing Strata
abn.pois <- glmer(Abundance ~ Tmean.day + Hmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(abn.pois)

# removing Hmean.day
abn.pois <- glmer(Abundance ~ Tmean.day + Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(abn.pois)

# removing Tmean.day
abn.pois <- glmer(Abundance ~ Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                    (1|SU), data = data_models, family = poisson)
vif(abn.pois)

# removing Hmean.day:Strata
abn.pois <- glmer(Abundance ~ Tsd.night + 
                    Tmean.day:Strata + Tsd.night:Strata + 
                    (1|SU), data = data_models, family = poisson)
vif(abn.pois) # all <3

summary(abn.pois)
deviance(abn.pois)/df.residual(abn.pois) # overdispersed

# considering negative binomial distribution
abn.nb <- glmer.nb(Abundance ~ Strata + Tmean.day + Hmean.day + Tsd.night + 
                     Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                     (1|SU), data = data_models)
vif(abn.nb)

# removing Hmean.day
abn.nb <- glmer.nb(Abundance ~ Strata + Tmean.day + Tsd.night + 
                     Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                     (1|SU), data = data_models)
vif(abn.nb)

# removing Strata
abn.nb <- glmer.nb(Abundance ~ Tmean.day + Tsd.night + 
                     Tmean.day:Strata + Tsd.night:Strata + Hmean.day:Strata +
                     (1|SU), data = data_models)
vif(abn.nb) # all <3

summary(abn.nb)

model.sel(abn.nb.null, abn.nb)
r.squaredGLMM(abn.nb)

saveRDS(ric.pois, here::here("output/output_model_ric.rds"))
saveRDS(abn.nb, here::here("output/output_model_abn.rds"))