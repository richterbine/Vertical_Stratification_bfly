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


#### usando as matrizes de particionamento de variancia Baselga 2013
dist <- bray.part(comm.spp)
dBC <- dist$bray
dBAL<- dist$bray.bal
dGRA<- dist$bray.gra

m_sd <- data.frame(Mean = c(mean(dBC, na.rm = T), mean(dBAL, na.rm = T), mean(dGRA, na.rm = T)),
                   SD = c(sd(dBC, na.rm = T), sd(dBAL, na.rm = T), sd(dGRA, na.rm = T)))

rownames(m_sd) <- c("dBC", "dBAL", "dGRA")

# preparation to run adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(comm.env, SU) # consider the "random effect" of Sampling unity

## model to evaluate the effect of covariates in the total dissimilarity
comm.m1 <- adonis2(dBC ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                  data = comm.env, permutations = perm, by = "marg")
comm.m1
# efeito do estrato sobre a variação na comunidade

### model to evauate the effect of covariates in the "turnover" component
comm.m2 <- adonis2(dBAL ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                  data = comm.env, permutations = perm, by="margin")
comm.m2
# efeito do estrato no turnover, canopy com menos turnover?

### model to evaluate the effect of covariates in the "nestedness" component
comm.m3 <- adonis2(dGRA ~ Tmean.day + Tsd.night + Hmean.day + Strata, 
                  data = comm.env, permutations = perm, by="marg")
comm.m3
# a variação na temperatura da noite afeta o quanto as espécies são aninhadas

# verify the homogeneity on the variance for each component
group <- as.factor(comm.env$Strata)

disp.dbc <- betadisper(dBC, group, type = "centr")
plot(betadisper(dBC, group, type = "centr"))
anova(betadisper(dBC, group, type = "centr"))
boxplot(betadisper(dBC, group, type = "centr"))


disp.dbal <- betadisper(dBAL, group, type = "centr")
plot(betadisper(dBAL, group, type = "centr"))
anova(betadisper(dBAL, group, type = "centr"))
boxplot(betadisper(dBAL, group, type = "centr"))


disp.dgra  <- betadisper(dGRA, group, type = "centr")
plot(betadisper(dGRA, group, type = "centr"))
anova(betadisper(dGRA, group, type = "centr"))
boxplot(betadisper(dGRA, group, type = "centr"))
# todos tem homogeneidade na variancia entre grupos

# PLOT BETADISPER GRAPHS --------------------------------------------------

### betadisper graphics
# for species data
rownames(comm.env)


data.disper <- data.frame(Component = c(rep("Total", nrow(comm.env)), 
                                        rep("Balance", nrow(comm.env)), 
                                        rep("Gradient", nrow(comm.env))),
                          Strata = rep(substr(rownames(comm.env), 1, 1), 3),
                          PCoA1 = c(disp.dbc$vectors[,"PCoA1"], disp.dbal$vectors[,"PCoA1"], 
                                    disp.dgra$vectors[,"PCoA1"]), 
                          PCoA2 = c(disp.dbc$vectors[,"PCoA2"], disp.dbal$vectors[,"PCoA2"],
                                    disp.dgra$vectors[,"PCoA2"]),
                          centroid1 = c(c(rep(disp.dbc$centroids[1, "PCoA1"], 38), rep(disp.dbc$centroids[2, "PCoA1"], 40)),
                                        c(rep(disp.dbal$centroids[1, "PCoA1"], 38), rep(disp.dbal$centroids[2, "PCoA1"], 40)),
                                        c(rep(disp.dgra$centroids[1, "PCoA1"], 38), rep(disp.dgra$centroids[2, "PCoA1"], 40))),
                          centroid2 = c(c(rep(disp.dbc$centroids[1, "PCoA2"], 38), rep(disp.dbc$centroids[2, "PCoA2"], 40)),
                                        c(rep(disp.dbal$centroids[1, "PCoA2"], 38), rep(disp.dbal$centroids[2, "PCoA2"], 40)),
                                        c(rep(disp.dgra$centroids[1, "PCoA2"], 38), rep(disp.dgra$centroids[2, "PCoA2"], 40))))

df <- subset(data.disper, Component == 'Total')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.bc <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  labs(tag = "a)")

df <- subset(data.disper, Component == 'Balance')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.bal <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  labs(tag = "b)")

df <- subset(data.disper, Component == 'Gradient')
find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- plyr::ddply(df, "Strata", find_hull)

plot.gra <- ggplot(data = df, aes(x = PCoA1, y = PCoA2, colour= Strata, fill = Strata)) +
  geom_point() + geom_point(aes(x= centroid1, y= centroid2), size= 3) +
  geom_segment(aes(x= centroid1, y= centroid2, xend= PCoA1, yend= PCoA2)) +
  geom_polygon(data = hulls, alpha = 0.2) + theme(legend.position = "none") +
  scale_color_manual(values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(values = c("firebrick3", "palegreen3")) +
  labs(tag = "c)")

# joint plots
plot.betadisper <- cowplot::plot_grid(plot.bc, plot.bal, plot.gra,
                   label_fontface = "plain", ncol = 3)


# relationship among environmental variables and groups composition
# only for data exploration
comm.env

label <- rownames(as.matrix(dist$bray.bal))
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

data <- data.frame(Total = as.numeric(dBC),
                   Balance = as.numeric(dBAL),
                   Gradient = as.numeric(dGRA),
                   pairs = nomes)

data$pairs <- as.factor(data$pairs)
levels(data$pairs)

dbc.st <- ggplot(data = subset(data, pairs != "UC"), aes(x= pairs, y= Total)) + 
  geom_boxplot(aes(colour = pairs, fill = pairs), alpha = .5) +
  #geom_smooth() + 
  labs(x = "Strata", y = expression(d[BC]), tag = "d)") +
  scale_color_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                     values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                     values = c("firebrick3", "palegreen3"))
dbc.st

dbal.st <- ggplot(data = subset(data, pairs != "UC"), aes(x= pairs, y= Balance)) + 
  geom_boxplot(aes(colour = pairs, fill = pairs), alpha = .5) +
  #geom_smooth() + 
  labs(x = "Strata", y = expression(d[BAL]), tag = "e)") +
  scale_color_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                     values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                    values = c("firebrick3", "palegreen3"))
dbal.st


Tsd.dist <- vegdist(comm.env$Tsd.night, method = "euclidean")
data$Tsd.night  <- as.numeric(Tsd.dist)
data %>% head

dgra.tsd <- ggplot(data = subset(data, pairs != "UC"), aes(x = Tsd.night, y = Gradient)) +
  geom_point(aes(colour = pairs), alpha = .5) +
  geom_smooth(color = "black") + labs(x = "Mean Humidity", y = expression(d[GRA]), tag = "f)") +
  scale_color_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                     values = c("firebrick3", "palegreen3")) +
  scale_fill_manual(name = "Site pairs", labels = c("Canopy - Canopy", "Understory - Understory"),
                    values = c("firebrick3", "palegreen3"))
dgra.tsd

plot.betadisper <- cowplot::plot_grid(plot.bc, plot.bal, plot.gra,
                                      dbc.st + theme(legend.position = "none"),
                                      dbal.st + theme(legend.position = "none"),
                                      dgra.tsd + theme(legend.position = "none"),
                                      label_fontface = "plain", ncol = 3)
plot.betadisper

cowplot::save_plot(here::here("output/Images/G_betadisper.png"), plot.betadisper,
                  base_height = 8, base_width = 12)
