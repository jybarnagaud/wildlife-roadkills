#--------------------------------------------------------------------#
#### This script computes the sampling effort covariate  ####

# for manuscript : 
# Spatial analysis of the wildlife roadkill risk at a regional scale. 
# Supplementary Materials
# anonymized until review process completed by the editor's request
#--------------------------------------------------------------------#

library(sf)
library(ade4)
library(ggplot2)
library(patchwork)
library(viridis)
library(factoextra)
library(adegraphics)

## data ------------------------------------------------------------------------

shape <- read_sf(dsn = "data/vector/maillage_5km_sampling_clean_12_05.shp")
shape1 <- as.data.frame(shape)

## retrieve sampling effort variables ------------------------------------------

# n.data = total number of records all species confounded
# n.dates = total number of unique dates with records all species confounded
# n.sp = species richness all taxa confounded
# n.obs = total number of observers

sampling <- data.frame(
  code.grid = shape1$Code, 
  name.grid = shape1$Nom,
  n.data = shape1$Nb_d_dn,
n.dates = shape1$Nb_d_dt,
n.sp = shape1$Nb_d_sp,
n.obs = shape1$Nb_d_bs)

sampling.pc <- sampling[,-c(1,2)]
colnames(sampling.pc) <- c("N. records", "N. dates", "N. species", "N. observers")

# retrieve geometry for later mapping

geom.shape <- st_geometry(shape)
centroid.shape <- st_centroid(geom.shape)
coords.shape <- st_coordinates(centroid.shape)

## compute the principal component analysis ------------------------------------

pc.sampeff <- dudi.pca(sampling.pc, nf = 2, scannf = F)

# % of explained variance

screeplot(pc.sampeff)
pc.sampeff$eig/sum(pc.sampeff$eig)

# correlation circle

s.corcircle(pc.sampeff$co)

# projection of points 

s.label(pc.sampeff$li)

## compute the sampling effort variable ----------------------------------------

pc1 <- -1*pc.sampeff$li[,1]
samp.eff <- sqrt(pc1+abs(min(pc1)))

## write results for analysis --------------------------------------------------

sampling2 <- cbind(coords.shape,sampling,samp.eff)
write.csv2(sampling2,"outputs/sampling_effort.csv")

## maps for Supplementary Materials -----------------------------------------

nbdd <- ggplot(shape)+
  geom_sf(aes(fill = log(Nb_d_dn+1)))+
  labs(fill = "Number of records \n (log-transformed)")+
  theme_classic()+
  scale_fill_viridis(option ="A")+ 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

nbdt <- ggplot(shape)+
  geom_sf(aes(fill = log(Nb_d_dt+1)))+
  labs(fill = "Number of dates \n (log-transformed)")+
  theme_classic()+
  scale_fill_viridis(option = "D")+ 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

nbsr <- ggplot(shape)+
  geom_sf(aes(fill = log(Nb_d_sp+1)))+
  labs(fill = "Number of species \n (log-transformed)")+
  theme_classic()+
  scale_fill_viridis(option = "E")+ 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

nbobs <- ggplot(shape)+
  geom_sf(aes(fill = log(Nb_d_bs+1)))+
  labs(fill = "Number of observers \n (log-transformed)")+
  theme_classic()+
  scale_fill_viridis(option = "G")+ 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

(nbdd + nbdt) / (nbsr + nbobs) 


ggsave("outputs/SM_2_maps.png",width=8, height = 8)

## PCA projections for SM ------------------------------------------------------

# normed scores
pc.sampeff$c1
