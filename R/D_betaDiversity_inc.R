# considering only the incidence of species
# BETA DIVERSITY ANALYSIS
library(betapart)
library(vegan)
library(ggplot2)
library(tidyverse)

# read data
comm.bfly <- readRDS(here::here("output/community_env_data.rds"))

# separate the community data (sites x Species) from environmental data (sites X env)
head(comm.bfly)
rownames(comm.bfly) <- paste(substr(comm.bfly$Strata, 1, 1), comm.bfly$Site, 
                             comm.bfly$SU, sep = "")
comm.bfly <- comm.bfly[-which(rowSums(comm.bfly[,4:41]) == 0),]

comm.spp <- (comm.bfly[,4:41])
comm.spp %>% head

comm.spp <- comm.spp[,-which(colSums(comm.spp) == 0)]

comm.env <- comm.bfly[,c(1:3, 43:50)]
comm.env %>% head()

inc.comm <- decostand(comm.spp, method = "pa")

bfly.dist <- beta.pair(inc.comm, index.family = "sor")
dT <- bfly.dist$beta.sor
dtur <- bfly.dist$beta.sim
dnes <- bfly.dist$beta.sne

M.SD <- data.frame(Mean = c(mean(dT, na.rm = T), mean(dtur, na.rm = T), mean(dnes, na.rm = T)),
           SD = c(sd(dT, na.rm = T), sd(dtur, na.rm = T), sd(dnes, na.rm = T)))

rownames(M.SD) <- c("dT", "dtur", "dnes")
M.SD

# preparation to run adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(comm.env, SU) # consider the "random effect" of Sampling unity

## model to evaluate the effect of covariates in the total dissimilarity
comm.m1.i <- adonis2(dT ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by = "marg")
comm.m1.i
# efeito do estrato sobre a variação na comunidade

### model to evauate the effect of covariates in the "turnover" component
comm.m2.i <- adonis2(dtur ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by="margin")
comm.m2.i
# efeito do estrato no turnover, canopy com menos turnover?

### model to evaluate the effect of covariates in the "nestedness" component
comm.m3.i <- adonis2(dnes ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by="marg")
comm.m3.i
# a variação na temperatura da noite afeta o quanto as espécies são aninhadas

# verify the homogeneity on the variance for each component
group <- as.factor(comm.env$Strata)

disp.dT <- betadisper(dT, group, type = "centr")
plot(betadisper(dT, group, type = "centr"))
anova(betadisper(dT, group, type = "centr"))
boxplot(betadisper(dT, group, type = "centr"))


disp.dtur <- betadisper(dtur, group, type = "centr")
plot(betadisper(dtur, group, type = "centr"))
anova(betadisper(dtur, group, type = "centr"))
boxplot(betadisper(dtur, group, type = "centr"))


disp.dnes  <- betadisper(dnes, group, type = "centr")
plot(betadisper(dnes, group, type = "centr"))
anova(betadisper(dnes, group, type = "centr"))
boxplot(betadisper(dnes, group, type = "centr"))

dist.tsd <- dist(comm.env$Tsd.night)

plot(dnes ~ dist.tsd)
abline(lm(dnes ~ dist.tsd))
