
# run IndVal and iNext analysis -------------------------------------------

# Required packages
#install.packages("indicspecies")
library(indicspecies)
library(tidyverse)

# Read data
comm.bfly <- readRDS(here::here("output/community_env_data.rds"))

# separate the community data (sites x Species) from environmental data (sites X env)
head(comm.bfly)
rownames(comm.bfly) <- paste(substr(comm.bfly$Strata, 1, 1), comm.bfly$Site, 
                             comm.bfly$SU, sep = "")
comm.bfly <- comm.bfly[-which(rowSums(comm.bfly[,4:41]) == 0),]
comm.bfly <- (comm.bfly[,4:41])

group <- as.factor(substr(rownames(comm.bfly), 1, 1))
group <- ifelse(group == "C", 1, 2)


# Indicator Value analysis
IndVal.bfly <- multipatt(comm.bfly, group, control = how(nperm = 999))
summary(IndVal.bfly, indvalcomp = TRUE)


# Performing the interpolation and extrapolation for diversity est --------

#install.packages("iNEXT")
library(iNEXT)

# Rarefection/Extrapolation curves
head(comm.bfly)

data.inext <- list(Canopy = colSums(comm.bfly[which(substr(rownames(comm.bfly), 1,1) == "C"),]),
                   Understory = colSums(comm.bfly[which(substr(rownames(comm.bfly), 1,1) == "U"),]))
data.inext %>% str

plot(sort(data.inext$Canopy, decreasing = T))
plot(sort(data.inext$Understory, decreasing = T))

DataInfo(data.inext)

out <- iNEXT(data.inext, q = c(0, 1, 2), nboot = 1000, datatype = "abundance")
out$iNextEst$Canopy$order <- as.factor(out$iNextEst$Canopy$order)
out$iNextEst$Understory$order <- as.factor(out$iNextEst$Understory$order)

levels(out$iNextEst$Canopy$order)[levels(out$iNextEst$Canopy$order) == "0"] <- "q = 0"
levels(out$iNextEst$Canopy$order)[levels(out$iNextEst$Canopy$order) == "1"] <- "q = 1"
levels(out$iNextEst$Canopy$order)[levels(out$iNextEst$Canopy$order) == "2"] <- "q = 2"

levels(out$iNextEst$Understory$order)[levels(out$iNextEst$Understory$order) == "0"] <- "q = 0"
levels(out$iNextEst$Understory$order)[levels(out$iNextEst$Understory$order) == "1"] <- "q = 1"
levels(out$iNextEst$Understory$order)[levels(out$iNextEst$Understory$order) == "2"] <- "q = 2"


# Sample-size-based R/E curves, separating by "site""
p.samp.size <- ggiNEXT(out, type=1, facet.var="order") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  theme(legend.position = "right")


# Coverage-based R/E curves, separating by "order"
p.cov.base <- ggiNEXT(out, type=3, facet.var="order")  +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  theme(legend.position = "right") 

RE_curves <- cowplot::plot_grid(p.samp.size + theme(legend.position = "none"),
                   p.cov.base + theme(legend.position = "none"))

cowplot::save_plot(here::here('output/Images/G_RE_curves.png'), RE_curves,
                   base_height = 6, base_width = 14)

# Diversity profile estimator
source(here::here("R/functions/ChaoHill_function.R"))
prof.can <- ChaoHill(data.inext$Canopy, datatype = "abundance")
prof.und <- ChaoHill(data.inext$Understory, datatype = "abundance")

rep(rownames(prof.can$EST), each = 2)

# joint results in a dataframe
data.profile <- data.frame(Values = c(as.vector(t(prof.can$EST)), as.vector(t(prof.und$EST))),
                           type = c(rep(rownames(prof.can$EST), each = ncol(prof.can$EST)), 
                                    rep(rownames(prof.und$EST), each = ncol(prof.und$EST))),
                           q = c(rep(colnames(prof.can$EST), 2), 
                                     rep(colnames(prof.und$EST), 2)),
                           LCI = c(as.vector(t(prof.can$LCI)), as.vector(t(prof.und$LCI))),
                           UCI = c(as.vector(t(prof.can$UCI)), as.vector(t(prof.und$UCI))),
                           Strata = c(rep("Canopy", 2*ncol(prof.can$EST)), rep("Understory", 2*ncol(prof.und$EST)))
                           )
head(data.profile)
data.profile$q <- substr(data.profile$q, 5,7)
data.profile$q <- as.numeric(data.profile$q)
data.profile$type <- ifelse(data.profile$type == "Observed", "Empirical", "Proposed")

profile <- ggplot(data = data.profile, aes(x = q, y = Values,  colour = Strata)) +
  geom_line(size = .8) + facet_wrap(~type) +
  geom_ribbon(aes(ymin= LCI, ymax= UCI,  fill = Strata),
                         alpha = .1, linetype = 0) +
  scale_x_continuous(breaks = c(0,1,2,3)) + scale_y_sqrt() +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  labs(x = "q order", y = "Hill numbers") 
profile

cowplot::save_plot(here::here("output/Images/G_profile_plot.png"), 
                   profile + theme(legend.position = "none"),
                   base_height = 4, base_width = 8)
