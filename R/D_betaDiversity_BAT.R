# BETA DIVERSITY ANALYSIS
library(betapart)
library(BAT)
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

# Beta diversity partitioning following the method proposed by Podani & Schmera (2011) and Carvalho et al. (2012)
beta_res <- BAT::beta(comm.spp, abund = T, runs = 499, func = "Soren")
Btot <- beta_res$Btotal
Brep <- beta_res$Brepl
Bric <- beta_res$Brich

m_sd <- data.frame(Mean = c(mean(Btot, na.rm = T), mean(Brep, na.rm = T), mean(Bric, na.rm = T)),
                   SD = c(sd(Btot, na.rm = T), sd(Brep, na.rm = T), sd(Bric, na.rm = T)))

rownames(m_sd) <- c("Btot", "Brep", "Bric")
round(m_sd, 3)

# preparation to run adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(comm.env, SU) # consider the "random effect" of Sampling unity

## model to evaluate the effect of covariates in the total dissimilarity
comm.m1 <- adonis2(Btot ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by = "marg")
comm.m1
# efeito do estrato sobre a variação na comunidade

### model to evauate the effect of covariates in the "turnover" component
comm.m2 <- adonis2(Brep ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by="margin")
comm.m2
# a slighter effect of strata but with a similar R2 value than mod1

### model to evaluate the effect of covariates in the "nestedness" component
comm.m3 <- adonis2(Bric ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                   data = comm.env, permutations = perm, by="marg")
comm.m3
# there aren't an effect of environemntal variables in abundance gain/loss among sites

# verify the homogeneity on the variance for each component
group <- as.factor(comm.env$Strata)

disp.Btot <- betadisper(Btot, group, type = "centr")
plot(betadisper(Btot, group, type = "centr"))
anova(betadisper(Btot, group, type = "centr"))
boxplot(betadisper(Btot, group, type = "centr"))


disp.Brep <- betadisper(Brep, group, type = "centr")
plot(betadisper(Brep, group, type = "centr"))
anova(betadisper(Brep, group, type = "centr"))
boxplot(betadisper(Brep, group, type = "centr"))


disp.Bric  <- betadisper(Bric, group, type = "centr")
plot(betadisper(Bric, group, type = "centr"))
anova(betadisper(Bric, group, type = "centr"))
boxplot(betadisper(Bric, group, type = "centr"))
# todos tem homogeneidade na variancia entre grupos

# PLOT BETADISPER GRAPHS --------------------------------------------------

### betadisper graphics
# for species data
rownames(comm.env)


data.disper <- data.frame(Component = c(rep("Total", nrow(comm.env)), 
                                        rep("Replacement", nrow(comm.env)), 
                                        rep("Richness", nrow(comm.env))),
                          Strata = rep(substr(rownames(comm.env), 1, 1), 3),
                          PCoA1 = c(disp.Btot$vectors[,"PCoA1"], disp.Brep$vectors[,"PCoA1"], 
                                    disp.Bric$vectors[,"PCoA1"]), 
                          PCoA2 = c(disp.Btot$vectors[,"PCoA2"], disp.Brep$vectors[,"PCoA2"],
                                    disp.Bric$vectors[,"PCoA2"]),
                          PCoA2 = c(disp.Btot$vectors[,"PCoA2"], disp.Brep$vectors[,"PCoA2"],
                                    disp.Bric$vectors[,"PCoA2"]),
                          centroid1 = c(c(rep(disp.Btot$centroids[1, "PCoA1"], 38), rep(disp.Btot$centroids[2, "PCoA1"], 40)),
                                        c(rep(disp.Brep$centroids[1, "PCoA1"], 38), rep(disp.Brep$centroids[2, "PCoA1"], 40)),
                                        c(rep(disp.Bric$centroids[1, "PCoA1"], 38), rep(disp.Bric$centroids[2, "PCoA1"], 40))),
                          centroid2 = c(c(rep(disp.Btot$centroids[1, "PCoA2"], 38), rep(disp.Btot$centroids[2, "PCoA2"], 40)),
                                        c(rep(disp.Brep$centroids[1, "PCoA2"], 38), rep(disp.Brep$centroids[2, "PCoA2"], 40)),
                                        c(rep(disp.Bric$centroids[1, "PCoA2"], 38), rep(disp.Bric$centroids[2, "PCoA2"], 40))))

df <- subset(data.disper, Component == 'Total')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.bc <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  annotate(geom = "text", x = -0.45, y = -0.45, size = 5, label = expression(beta[tot])) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  labs(tag = "a)")
plot.bc

# beta replacement
df <- subset(data.disper, Component == 'Replacement')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.bal <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  annotate(geom = "text", x = -0.38, y = -0.38, size = 5, label = expression(beta[rep])) +
  labs(tag = "b)")
plot.bal

# beta abundance difference
df <- subset(data.disper, Component == 'Richness')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.gra <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  annotate(geom="text", x= -.42, y= -.16, label= expression(beta[ric]), size = 5) +
  labs(tag = "c)")
plot.gra

# relationship among environmental variables and groups composition
# only for data exploration
comm.env

label <- rownames(as.matrix(beta_res$Brepl))
label.col <- label[-78]
label.row <- label[-1]

nomes <- as.vector(NA)

for (i in 1:length(label.col)) {
  if(i == 1) {
    nomes <- paste(substr(label.row, 1, 1), substr(label.col[i], 1, 1), sep = "")
  } else{
    nomes <- c(nomes, paste(substr(label.row[-c(1:i-1)], 1, 1), substr(label.col[i], 1,1), sep = ""))
  }
}

data <- data.frame(Total = as.numeric(Btot),
                   Replacement = as.numeric(Brep),
                   Richness = as.numeric(Bric),
                   pairs = nomes)

data$pairs <- as.factor(data$pairs)
levels(data$pairs)

comm.m1
comm.m2

Btot.st <- ggplot(data = subset(data, pairs != "UC"), aes(x= pairs, y= Total)) + 
  geom_boxplot(aes(colour = pairs, fill = pairs), alpha = .5) +
  #geom_smooth() + 
  labs(x = "Strata", y = expression(beta[tot]), tag = "d)") +
  scale_color_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                     values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                    values = c("firebrick3", "palegreen3"))
Btot.st


plot.betadisper <- cowplot::plot_grid(plot.bc, plot.bal, plot.gra,
                                      Btot.st + theme(legend.position = "none"),
                                      label_fontface = "plain", ncol = 2)
plot.betadisper

cowplot::save_plot(here::here("output/Images/Fig4_betadisper.png"), plot.betadisper,
                   base_height = 6, base_width = 8)
