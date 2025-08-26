#-----------------------------------------------------------#
#### This script computes the sampling effort covariate  ####
#-----------------------------------------------------------#

library(sf)
library(ade4)

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


