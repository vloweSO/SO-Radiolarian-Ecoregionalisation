# Ecoregionalisation of the Southern Ocean using Radiolarians
# Adapted from 'Numerical Ecology with R' Borcard, Gillet and Legendre 2018
# Vikki Lowe
# 8 March 2021
#==========================================================================================================================================
# LOAD ALL LIBRARIES, FUNCTIONS & DATA
#==========================================================================================================================================
library(dplyr)
library(ggplot2)
library(ggrepel)
library(labdsv)
library(mvpart)
library(MVPARTwrap)
library(vegan)

### Load FUnctions
source("drawmap.R")
source("triplot.rda.R")

### Load Data - Data should be changed 
spe <- read.csv("SOspe.csv")
spa <- read.csv("SOspa.csv")
factors <- read.csv("SOfactors.csv")
env <- read.csv("SOenv.csv")
spe.xy <- spa

row.names(spe) <- factors$Site
row.names(factors) <- factors$Site
row.names(spa) <- factors$Site
row.names(env) <- factors$Site
row.names(spe.xy) <- factors$Site

#==========================================================================================================================================
# DATA TRANFORMATIONS
#==========================================================================================================================================
spe.hel <- decostand(spe, "hellinger")
spe.norm <- decostand(spe, "normalize")

#==========================================================================================================================================
# NMDS
#==========================================================================================================================================

datanMDS <- metaMDS(spe.hel, k=2, trymax = 999, distance = "euclidean" , wascores = TRUE)
Data.scores = as.data.frame(scores(datanMDS))
Data.scores = mutate(Data.scores, data = seq(n()))
factors.scores = mutate(factors, data = seq(n()))
Data.scores = right_join(factors.scores, Data.scores)
datanMDS$stress <- round(datanMDS$stress, 4)

Data.scores$MRTClust <- as.factor(Data.scores$MRTClust)
plot = ggplot() +
  geom_point(data=Data.scores, aes(x=NMDS1, y=NMDS2, shape = Zone, colour = Sector), size = 3) +
  #geom_text_repel(data=Data.scores, aes(x=NMDS1, y=NMDS2, label=Site), size = 3) + 
  scale_shape_manual(values = c(18, 15, 16, 17, 20, 21, 19)) +
  scale_color_manual(values = c("green3", "darkorange", "purple", "red", "black", "grey", "blue", "yellow")) +
  #scale_color_manual(values = c("red", "darkgreen", "turquoise", "purple", "brown", "orange", "blue", "green", "grey81", "grey16", "yellow", "deeppink")) +
  annotate("text", x = -0.4, y = 0.6, label = paste("Stress = ", datanMDS$stress)) +
  ggtitle("SO") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot(plot)

#==========================================================================================================================================
# RDA
#==========================================================================================================================================

spe.rda <- rda(spe.hel ~ ., env, na = na.omit)
summary(spe.rda)
(R2 <- RsquareAdj(spe.rda)$r.squared)
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

PARTspe.rda <- rda(spe.hel ~ T_200 + T_300 + T_400 + T_500 + N_200 + Depth_27.9 + O_10 + Den_500 + N_300 + T_100 + P_200 + T_10 + N_100 + Depth_27.3 + P_100 + Depth_27.6 + Den_200 + N_400 + Si_500 + P_300 + Den_10 + P_10 + N_10 + Si_300 + P_400 + Si_200 + Den_100 + N_500, env, na = na.omit)
(R2 <- RsquareAdj(PARTspe.rda)$r.squared)
(R2adj <- RsquareAdj(PARTspe.rda)$adj.r.squared)
summary(PARTspe.rda)

### Species at least 0.2 in ordination plane of axes 1 & 2
spe.good <- goodness(spe.rda)
sel.sp <- which(spe.good[, 2] >= 0.2)

PARTspe.good <- goodness(PARTspe.rda)
PARTsel.sp <- which(PARTspe.good[, 2] >= 0.2)

### Triplots for full and partial RDA
dev.new(
  title = "RDA plot with triplot.rda - ALL SORAD",
  width = 20,
  height = 16,
  noRStudioGD = TRUE)
#par(mfrow = c(1, 2))
triplot.rda(spe.rda,
            site.sc = "lc", 
            scaling = 1, 
            cex.char2 = 0.7, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05, 
            select.spe = sel.sp
)

dev.new(
  title = "PARTRDA plot with triplot.rda - SO PART",
  width = 20,
  height = 16,
  noRStudioGD = TRUE)
#par(mfrow = c(1, 2))
triplot.rda(PARTspe.rda,
            site.sc = "lc", 
            scaling = 1, 
            cex.char2 = 0.7, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05, 
            select.spe = PARTsel.sp
)

#==========================================================================================================================================
#GLOBAL TEST FOR RDA & TEST CANONICAL AXES
#==========================================================================================================================================

anova(spe.rda, permutations = how(nperm = 999))
anova(spe.rda, by = "axis", permutations = how(nperm = 999))

#==========================================================================================================================================
# MRT : CONSTRAINED CLUSTERING
#==========================================================================================================================================

#ALL
dev.new(
  title = "Multivariate regression tree - all explanatory variables - ALL",
  width = 14,
  height = 7,
  noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
spe.ch.mvpart <-
  mvpart(
    data.matrix(spe.norm) ~ .,
    env,
    margin = 0.08,
    cp = 0,
    xv = "pick",
    xval = nrow(spe),
    xvmult = 100
  )

summary(spe.ch.mvpart)
printcp(spe.ch.mvpart)

dev.new(
  title = "Residuals of MRT",
  width = 10,
  height = 6,
  noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
hist(residuals(spe.ch.mvpart), col = "bisque")
plot(predict(spe.ch.mvpart, type = "matrix"),
     residuals(spe.ch.mvpart),
     main = "Residuals vs Predicted")
abline(h = 0, lty = 3, col = "grey")


spe.ch.mvpart.wrap <-
  MRT(spe.ch.mvpart, percent = 10, species = colnames(spe))
summary(spe.ch.mvpart.wrap)
spe.ch.MRT.indval <- indval(spe.norm, spe.ch.mvpart$where)
pval.adj3 <- p.adjust(spe.ch.MRT.indval$pval)
spe.ch.MRT.indval$maxcls[which(pval.adj3 <= 0.05)]
spe.ch.MRT.indval$indcls[which(pval.adj3 <= 0.05)]
spech.mvpart.g <- factor(spe.ch.mvpart$where)
levels(spech.mvpart.g) <- 1:length(levels(spech.mvpart.g))
table(spech.mvpart.g, spech.UPGMA.g)

dev.new(title = "MRT clusters",
        width = 9,
        noRStudioGD = TRUE)
drawmap(xy = spa,
        clusters = spech.mvpart.g,
        main = "MRT clusters")