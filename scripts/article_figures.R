#------------------------------------------------------------------------------#
##### Figures for [XXX et al.] #####
# corresponding author : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
#------------------------------------------------------------------------------#
                    
## libraries -------------------------------------------------------------------

library(sf)
library(inlabru)
library(INLA)
library(fmesher)
library(tmap)
library(raster)
library(ade4)
library(factoextra)
library(adegraphics)
library(patchwork)
library(ggpubr)
library(tidyr)
library(forestplot)

# additional functions
source("functions/chargeDON_observations.R")
source("functions/inlabookfunctions.R")

## data ------------------------------------------------------------------------

# covariates, point patterns and utilities for the model

load("data/OUT_chargeDON_covariables_mesh.rda")

# general contour for France

france.contour <- st_read(dsn = "data/France_contour.geojson")

# contour of the study area

empr <- read_sf("data/domaine/emprise_maillage_1km_g2.shp")
empr2 <- st_geometry(empr)
class(empr2)
empr.sf <- st_as_sf(empr2)

# model fits

fit <- read.csv2("outputs/ajustement.csv", sep = " ", dec = ".")  # note : this is an output 

# model parameters

load("data/estim_effects.rda")

## SM 1: map of covariates =====================================================

# raster names and classes

lc.list <- data.frame(
  file = c(
    "urbain_landcover_1km.tif",
    "foret_landcover_1km.tif",
    "culture_landcover_1km.tif",
    "densite_hydro_1km.tif",
    "densite_haie_1km.tif"
  )
)

lc.list$type <- "lc"

raster.dist <- data.frame(
  file =
    c("raster_hydro_dist_50m.tif",
      "raster_veg_dist_1km.tif"),
  type = rep("dist", 2)
)

raster.road <- data.frame(file =
                            c("vitesse_1km.tif" ,
                              "trafic_1km.tif"),
                          type = "road")

raster.effort <- c("sampling_axis1_1km.tif", "sampling")

raster.data <- rbind(lc.list, raster.dist, raster.road)

raster.data$pal <-
  c("Greys",
    "Greens",
    "YlOrBr",
    "Blues",
    "BrBG",
    "Blues",
    "Greens",
    "plasma",
    "viridis")

raster.data$tit <-
  c(
    "Urban cover, ha",
    "Forest cover, ha",
    "Agriculture cover, ha",
    "Wetlands cover, ha",
    "Hedgerows, ha",
    "distance to water, m",
    "distance to vegetation, m",
    "speed, km/h",
    "traffic, cars/day"
  )

rownames(raster.data) <- c("urb","woo","agr","wet","hedg","diw","div","spe","traf")

tab.1 <- data.frame()

for(i in  1:nrow(raster.data)){
  
# load raster
  
ras <- raster(paste("data/raster_cov/", raster.data[i, "file"], sep = ""))

# for land covers: convert to hectares per km²

if (raster.data[i, "type"] == "lc") {
  ras2 <- ras / 10000
} else
  ras2 <- ras

# for speed and traffic : if null, set to NA

if(raster.data[i,"type"] == "road"){
  ras2 <- ras
  ras2[ras2 == 0] <- NA
}

# plot raster

raster.map <-
  tm_shape(ras2, bbox = st_bbox(c(
    xmin = -5.3,
    xmax = 0,
    ymax = 48.8,
    ymin = 46.8
  ), crs = st_crs(4326))) +
  tm_raster(
    palette = raster.data[i, "pal"],
    title = raster.data[i, "tit"],
    colorNA = "white",
    style = "pretty"
  ) +
  tm_layout(legend.position = c("left", "bottom")) +
  tm_shape(france.contour) +
  tm_borders(col = "black", alpha = 0.2)

assign(rownames(raster.data)[i], raster.map)

tmap_save(raster.map, paste("outputs/SM1_",rownames(raster.data)[i],".png",sep=""), dpi = 600)

# data for table 1 (NA converted to 0 : absence of given landscape feature)

z.val <- getValues(ras)
z.val.0 <- z.val
z.val.0[is.na(z.val.0)] <- 0

tab.1 <-
  rbind(tab.1,
        c(mean(z.val.0,na.rm=T), sd(z.val.0,na.rm=T), min(z.val.0,na.rm=T), max(z.val.0,na.rm=T),
        sum(z.val.0 > 0,na.rm=T),
        length(z.val.0)))

}

## TABLE 1 - description of raw variables --------------------------------------

tab.1b <- cbind(raster.data$tit,round(tab.1,2))
colnames(tab.1b) <- c("variable","mean","sd","min","max","nb>0","nb")

# note : computations are done in the previous section (SM1)

# add sampling data

samp.gis <- st_read("data/sampling/maillage_5km_sampling_clean_12_05.shp")

samp.data <- data.frame(
n.data = samp.gis$Nb_d_dn,
n.sp = samp.gis$Nb_d_sp,
n.obs = samp.gis$Nb_d_bs,
n.dates = samp.gis$Nb_d_dt
)
rownames(samp.data) <- samp.gis$id_plyg

mean.samp <- apply(samp.data,2,mean)
sd.samp <- apply(samp.data,2,sd)
min.samp <- apply(samp.data,2,min)
max.samp <- apply(samp.data,2,max)
nb0.samp <- apply(samp.data,2,function(x){sum(x>0)})
length.samp <- apply(samp.data,2,length)
vari.names <- c("sampling pressure","sampling diversity","number of observers","sampling frequency")

samp.df <- data.frame(vari.names,mean.samp,sd.samp,min.samp,max.samp,nb0.samp,length.samp)
colnames(samp.df) <- colnames(tab.1b)

# compile all 

tab.1c <- rbind(tab.1b,samp.df)
tab.1c$prop.notnull <- 100 * (tab.1c$'nb>0' / tab.1c$nb)
# write.csv2(tab.1c,"outputs/Table_1.csv")

## SM1 : land cover PCA --------------------------------------------------------

rownames(pc$co) <- c("forest", "cultures","wetlands","hedgerows","urban")

# axes 1 vs 2
p12 <- fviz_pca_biplot(pc, axes = c(1,2),geom = "point",
                col.ind = "gray80", fill.ind = "gray90",alpha.ind = 0.5,
                col.var = "#440154",title = "",)
p12b <- ggpar(p12, xlab = "PC1", ylab = 'PC2') 

# axes 2 vs 3
p23 <- fviz_pca_biplot(pc, axes = c(2,3),geom = "point",
                       col.ind = "gray80", fill.ind = "gray90",alpha.ind = 0.5,
                       col.var = "#440154",title = "",)
p23b <- ggpar(p23, xlab = "PC2", ylab = 'PC3') 

# axes 3 vs 4
p34 <- fviz_pca_biplot(pc, axes = c(3,4),geom = "point",
                       col.ind = "gray80", fill.ind = "gray90",alpha.ind = 0.5,
                       col.var = "#440154",title = "",)
p34b <- ggpar(p23, xlab = "PC3", ylab = 'PC4') 

(p12b + p23b) / (p34b + plot_spacer()) + plot_annotation(tag_levels = "a")

ggsave("outputs/SM1_lcPCA.png", height = 10, width = 10)

## SM1 : maps of sampling variables --------------------------------------------

plot(samp.gis["Nb_d_sp"])


## TABLE 2 - model fit ---------------------------------------------------------

# names to english 

rownames(fit) <-
  c(
    "wild boar",
    "red fox",
    "pine marten",
    "european badger",
    "roe deer",
    "house marten",
    "least weasel",
    "stoat"
  )

# retrieve some features

rownames(fit)[apply(fit,2,which.max)]

## FIGURE 3 - example of model fit ---------------------------------------------

species <- "Chevreuil europeen"
load("data/AUC_TSS_PRED_MORT.rda")

auc <- AUCtot[[species]]
tss <- TSStot[[species]]
pred <- PREDtot[[species]]
n.kill <- NobsMOR_res0[[species]]

pred.stat <- apply(pred,2,quantile,probs=c(0.025,0.5,0.975))

pred.stat1 <- data.frame(n.kill,t(pred.stat))
colnames(pred.stat1) <- c("n.kill","q1","q2","q3")

fit <- data.frame(auc,tss)

auc.chev <- ggplot(fit) + 
  aes(x = auc)+
  geom_histogram(bins = 10, fill = "gray90",col="gray30")+
  theme_classic()+
  labs(x = "AUC",y = "Frequency")
  
tss.chev <- ggplot(fit) + 
  aes(x = tss)+
  geom_histogram(bins = 10, fill = "gray90",col="gray30")+
  theme_classic()+
  labs(x = "TSS",y = "Frequency")

pred.error <- ggplot(pred.stat1)+
  aes(x = n.kill, y = q2)+
  geom_point()+
  geom_errorbar(aes(ymin =q1, ymax = q3,width = 0.05))+
  geom_abline(slope=1, intercept = 0,linetype = "dashed", col = "steelblue",linewidth = 1)+ 
  geom_smooth(method='lm', col = "goldenrod",linewidth =1 , se = F, linetype = "dashed")+
  theme_classic()+
  labs(x = "Number of observed roadkills",y = "Number of predicted roadkills")

(auc.chev + tss.chev) / (pred.error + plot_spacer()) + plot_annotation(tag_levels = "a")

# ggsave("outputs/fit_chevreuil.png", width = 10, height = 10)

# NOTE : still need to add effect of spatial scale

## FIGURE 4 - predicted maps for exposure layer - roe deer----------------------

# roe deer taken for display. shows raw data, predictions and spatial field

species <- "Chevreuil europeen"

raw.dat <- chargeDON(species,2015,2020,CRSdon=CRS(proj4string(foretSG)))

dat.liv <- raw.dat[[1]]
dat.liv.sf <- st_as_sf(dat.liv)
dat.kill <- raw.dat[[2]]
N.liv <- raw.dat[[3]]
N.kill <- raw.dat[[4]]

# prediction of the exposure layer : predicted distribution of alive animals

load(paste0("data/fitVIV_",species,".rda"))

pred.liv <- log(as.vector(predVIV[2,]))

projgrid <- inla.mesh.projector(mesh, coordinates(pxlVIV))

Wmean <- inla.mesh.project(projgrid, ppVIV$summary.random$i$mean) # extract random spatial field

v1 <- ggplot()+
  gg(pxlVIV, aes(fill=pred.liv))+
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "covariates + spatial field")+
  scale_fill_viridis_c(option = "viridis", name = "intensity")+
  theme_classic()+
  theme(axis.line=element_blank())

v2 <- ggplot() +
  gg(pxlVIV, aes(fill = Wmean)) +
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "spatial field only")+
  scale_fill_viridis_c(option = "inferno", name = "intensity")+
  theme_classic()+
  theme(axis.line=element_blank())

v3 <- ggplot() +
  gg(pxlVIV, aes(fill = pred.liv - Wmean)) +
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "covariates only")+
  scale_fill_viridis_c(option = "mako", name = "intensity")+
  theme_classic()+
  theme(axis.line=element_blank())

crs.inla <- crs(pxlVIV)
empr.sf2 <- empr.sf
empr.sf2 <- st_transform(empr.sf2, crs = st_crs(pxlVIV))

v4 <- ggplot(dat.liv.sf) +
  geom_sf(size = 0.5,col="black",alpha = 0.1)+
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "records")+
  theme_classic()+
  theme(axis.line=element_blank())

(v1 + v2) / (v3 + v4) + 
  plot_annotation(tag_levels = "a")

# ggsave("outputs/figure4_chevreuil.png")


## FIGURE 5 - predicted maps for risk layer - roe deer--------------------------

load(paste0("data/pred_",species,".rda"))

# prediction of intensity of killed animals given exposure

killed <- apply(matrix(unlist(lapply(predMORT, function(x)
  as.vector(x))), ncol = 100), 1, mean)

k1 <-
  ggplot() +
  gg(pxlMOR, aes(fill = log(killed)))+
  geom_sf(data = empr.sf2,alpha = 0,col="gray70") +
  labs(x = "longitude", y = "latitude", title = "covariates + spatial field")+
  scale_fill_viridis_c(option = "viridis", name = "intensity")+
  theme_classic()+
  theme(axis.line=element_blank())

# spatial field = residuals of intensity of killed animals

load(paste0("data/modINLA_",species,"1.rda"))

projgrid <- inla.mesh.projector(meshROUTE, coordinates(pxlMOR))

Wmean <- inla.mesh.project(projgrid, modINLA$summary.random$i$mean)

k2 <- 
  ggplot()+
  gg(pxlMOR, aes(fill=Wmean))+
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "spatial field")+
  scale_fill_viridis_c(option = "inferno", name = "intensity")+
  theme_classic()+
  theme(axis.line=element_blank())

# raw point pattern

k3 <- 
  ggplot()+
  gg(dat.kill,size=0.5,col="black")+
  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
  labs(x = "longitude", y = "latitude", title = "records")+
  theme_classic()+
  theme(axis.line=element_blank())


(k1 + k2) / (k3 + plot_spacer()) + 
  plot_annotation(tag_levels = "a")

#ggsave("outputs/figure5_chevreuil.png")

## SM2 - as figure 5, all species ----------------------------------------------

# [not run at the moment, missing data]

# for (i in 1:length(listeSP)){
#species <- listeSP[i]

#load(paste0("data/pred_",species,".rda"))

# prediction of intensity of killed animals given exposure

#killed <- apply(matrix(unlist(lapply(predMORT, function(x)
 # as.vector(x))), ncol = 100), 1, mean)

#k1 <-
 # ggplot() +
  #gg(pxlMOR, aes(fill = log(killed)))+
  #geom_sf(data = empr.sf2,alpha = 0,col="gray70") +
  #labs(x = "longitude", y = "latitude", title = "covariates + spatial field")+
  #scale_fill_viridis_c(option = "viridis", name = "intensity")+
  #theme_classic()+
  #theme(axis.line=element_blank())

# spatial field = residuals of intensity of killed animals

#load(paste0("data/modINLA_",species,"1.rda"))

#projgrid <- inla.mesh.projector(meshROUTE, coordinates(pxlMOR))

#Wmean <- inla.mesh.project(projgrid, modINLA$summary.random$i$mean)

#k2 <- 
#  ggplot()+
# gg(pxlMOR, aes(fill=Wmean))+
#  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
#  labs(x = "longitude", y = "latitude", title = "spatial field")+
#  scale_fill_viridis_c(option = "inferno", name = "intensity")+
#  theme_classic()+
#  theme(axis.line=element_blank())

# raw point pattern

#k3 <- 
#  ggplot()+
#  gg(dat.kill,size=0.5,col="black")+
#  geom_sf(data = empr.sf2,alpha = 0,col="gray70")+
#  labs(x = "longitude", y = "latitude", title = "records")+
#  theme_classic()+
#  theme(axis.line=element_blank())

#plot.i <-(k1 + k2) / (k3 + plot_spacer()) + 
#  plot_annotation(tag_levels = "a")

#ggsave(paste("outputs/SM2_",species,".png",sep=""),plot = plot.i)
#}

## FIGURE 6 - parameter values and confidence intervals ------------------------

# the forest plots display parameter values and 95% credibility intervals
# correspond to the average over 99 draws

# risk layer - all variables

list.kill.variables <- list('(a) Exposure' = Efaune,'(b) Distance to water' = Edisthydro,
                            '(c) Distance to vegetation' = Edistveg,
                             '(d) Speed' = Evitesse,
                            '(e) Traffic' = Etrafic
                            )

param.faune <- data.frame()

for(i in 1:length(list.kill.variables)){
 tp<-
  as.data.frame(t(sapply(
    list.kill.variables[[i]],
    FUN = function(x) {
      apply(x, 2, "mean", na.rm = T)
    }
  )))
 tp$vari <- names(list.kill.variables)[i]
 param.faune <- rbind(param.faune,tp)
}

param.faune$species <- c(
  "wild boar",
  "red fox",
  "pine marten",
  "stoat",
  "house marten",
  "roe deer",
  "european badger",
  "least weasel"
)

colnames(param.faune) <- c( "lower", "median", "upper","variable","species")
param.faune$col.sig <- "gray30"
param.faune[which(param.faune$lower > 0 &
                    param.faune$upper > 0), "col.sig"] <- "steelblue"
param.faune[which(param.faune$lower < 0 &
                    param.faune$upper < 0), "col.sig"] <- "darkred"

fp.faune <- ggplot(param.faune) +
  aes(y = species, x = median, col = col.sig) +
  geom_point(size = 2) +
  geom_linerange(aes(xmin = lower, xmax = upper),linewidth = 0.9) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             col = "gray70") +
  theme_classic() +
  theme(strip.text = element_text(face="bold", size=10))+
  scale_colour_identity()+
  labs(x = "Estimate - 95% CI", y = "Species")+
  facet_wrap(~variable,axes = "all")+ 
  theme(strip.background=element_blank())
fp.faune

#ggsave("outputs/Figure6.png")


## FIGURE XX - sampling effort -------------------------------------------------


## recomputing sampling effort from raw data -----------------------------------

## open data files

mmssp <-
  read.csv2(
    "D:/encadrement/Vandroux_2022/data/data/data_covariable/sampling/exportFB_Mam.csv"
  )

for(i in c(22,29,35,56)){
nm <- paste("avissp",i,sep="") 
fil <- paste("D:/encadrement/Vandroux_2022/data/data/data_covariable/sampling/exportFB_Avi",i,
             "_03052022.txt",sep="")
  assign(nm,
  read.table(fil,
    header = T,
    sep = "\t"
  )
  )
}

avissp <- rbind(avissp22, avissp29,avissp35,avissp56)

## recompute sampling effort for mammals

# species richness, number of observers, number of unique dates (birds and mammals confounded)

comp.effort <- function(x,y){

  tp <-
    unique(x[, c("GRID_NAME_maille",y)])

  tp.agr <-
    aggregate(
      tp[,y],
      by = list(tp$GRID_NAME_maille),
      FUN = "length")
}

vari <- c("ID_SPECIES_Biolovision","DATE","UNIVERSAL_ID_OBSERVER")
names(vari) <- c("sr","dat","obs")
dataset <- c("mmssp","avissp")

for(i in 1:3){
  dftp <- data.frame()
  for(j in 1:2){
    dtp <- get(dataset[j])
      atp <- comp.effort(dtp,vari[i])
      colnames(atp) <- c("GRID_NAME_maille","z")
    dftp <- rbind(dftp,atp)
  }
  sftp <- aggregate(dftp$z,by = list(dftp$GRID_NAME_maille),FUN="sum")
  colnames(sftp) <- c("GRID_NAME_maille",names(vari)[i])
  assign(paste("samp.new",names(vari)[i],sep="_"),sftp)
  }

samp.new.1 <- merge(samp.new_dat,samp.new_obs,by = "GRID_NAME_maille")
samp.new.n <- merge(samp.new.1,samp.new_sr,by = "GRID_NAME_maille")

# number of data, birds and mammals confounded

n.obs.mm <- aggregate(mmssp$UNIVERSAL_ID_OBSERVATION,by = list(mmssp$GRID_NAME_maille), FUN = "length")
n.obs.avi <- aggregate(avissp$UNIVERSAL_ID_OBSERVATION,by = list(avissp$GRID_NAME_maille), FUN = "length")
n.obs.tot <- merge(n.obs.mm,n.obs.avi,by = "Group.1")
n.obs.tot$n.obs <- apply(n.obs.tot[,2:3],1,"sum")
samp.new.all <- merge(n.obs.tot[,c(1,4)],samp.new.n,by.x = "Group.1",by.y = "GRID_NAME_maille")
colnames(samp.new.all)[1] <- "GRID_NAME_maille"

# compute a PCA axis

samp.new.pc <- samp.new.all[,-1]
rownames(samp.new.pc) <- samp.new.all[,1]
pc.new <- dudi.pca(samp.new.pc,scannf = F, nf = 2)
samp.new.all2 <- merge(samp.new.all,pc.new$li,by.x = "GRID_NAME_maille",by.y = 0)

# add coordinates

xy.samp <-
  aggregate(avissp[, c("COORD_LAT_.WGS84.", "COORD_LON_.WGS84.")],
            by = list(avissp$"GRID_NAME_maille"),
            FUN = "mean")
colnames(xy.samp) <-
  c("GRID_NAME_maille", "COORD_LAT_.WGS84.", "COORD_LON_.WGS84.")


samp.new.xy <-
  merge(samp.new.all2, xy.samp, by = "GRID_NAME_maille", all = T)

colnames(samp.new.xy) <-
  c("GRID_NAME_maille", "nb.data","nb.dates","nb.obs","sr" ,"Axis1","Axis2","lat.wgs84", "lon.wgs84")

# pass to sf

samp.sf <- st_as_sf(x = samp.new.xy, 
                    coords = c("lon.wgs84","lat.wgs84"),
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

plot(samp.sf)

# st_write(samp.sf, "outputs/new_sampling_raw.shp")

# save also the df
 
#write.csv2(samp.new.xy,"outputs/new_sampling_variables.csv")


