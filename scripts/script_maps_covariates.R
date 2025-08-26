#--------------------------------------------------------------------#
#### This script computes the maps and PCA for environmental covariates

# for manuscript : 
# Spatial analysis of the wildlife roadkill risk at a regional scale. 
# anonymized until review process completed by the editor's request
# Supplementary Materials 
# anonymized until review process completed by the editor's request
#--------------------------------------------------------------------#

library(raster)
library(terra)

# names of variables for legends

all.maps <- list.files("data/raster")
names(all.maps) <- c(
  "Agriculture cover, ha",
  "Hedgerows, m.ha-1",
  "Wetlands cover, ha",
  "Forest cover, ha",
  "Distance to water, m",
  "Distance to vegetation, m",
  "Sampling axis",
  "Traffic, cars.day-1",
  "Urban cover, ha",
  "Speed, km.h-1"
)



# loop over all variables

for (i in 1:length(all.maps)) {
  
  # open raster
  r1 <- raster(paste("data/raster/", all.maps[i], sep = ""))
  
  # for column names with superscripts in legends
  parts <- strsplit(names(all.maps)[i], "-")[[1]]
  base <- parts[1]
  expo <- parts[2]
    png(paste("outputs/map_", names(all.maps)[i], ".png", sep = ""))
  if (!is.na(expo)) {
    plot(r1, main = bquote(.(base)^.(paste0("-", expo))))
  } else
    (plot(r1, main = names(all.maps)[i]))
  dev.off()
}

## PCA on environmental variables ----------------------------------------------

# Note : this data table is generated from the maps above in the main script. To
# replicate it, rasters have to be clipped to the region of interest (Brittany 
# only, excludes the SE part of the raster). 

env.dat <- read.table("outputs/environment_data.txt",
                      header = T,
                      sep = " ")

pc.env <- dudi.pca(env.dat,scannf = F, nf = 4)

pc.env$eig/sum(pc.env$eig)
write.csv2(pc.env$c1,"outputs/SM_environment_loadings.csv")
s