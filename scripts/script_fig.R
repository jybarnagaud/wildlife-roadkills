library(sf)
library(fmesher)
library(inlabru)
library(INLA)
library(ggplot2)
library(igraph)
library(dplyr)
library(raster)
library(tidygraph)
library(tidyverse)
library(ade4)
library(fields)
library(pROC)

source("functions/inlabookfunctions.R")
#source("functions/chargeDON_covariables_mesh_modifJYB31072024.R")
load("data/OUT_chargeDON_covariables_mesh.rda")
source("functions/chargeDON_observations.R")

## Figure 1 = mesh + données ---------------------------------------------------

listeSP=intersect(unique(roadkill_sp$spNEW),unique(datVIV$nom_vern))
listeSP=listeSP[c(1:5,7:9)]
listeSP

for (espece in listeSP){
  DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
  datviv = DON[[1]]
  datcoll = DON[[2]]
  NobsVIV = DON[[3]]
  NobsMOR = DON[[4]]
  
  png(file=paste0("don_",espece,".png"),width=300,height=200,res=300,units="mm")
  par(mar=rep(1,4),mfrow=c(1,1))
  
  class(mesh) <- "fm_mesh_2d"
  plot(mesh, main = "", asp = 1)
  plot(routesSL,add=T,col=2,border=2,lwd=2)
  plot(domaineSP,add=T,border=1,lwd=2)
  plot(datviv,add=T,col=3,pch=16,cex=0.5)
  plot(datcoll,add=T,col=4,pch=16,cex=0.5)
  dev.off()
}

ggplot() +
  gg(mesh)

#fig3
espece="Sanglier"
load("data/AUC_TSS_PRED_MORT.rda")
AUC=AUCtot[[espece]]
TSS=TSStot[[espece]]
PRED=PREDtot[[espece]]
NobsMOR=NobsMOR_res0[[espece]]
PRED.stat=apply(PRED,2,quantile,probs=c(0.025,0.5,0.975))
par(mfrow=c(1,3))
hist(AUC)
hist(TSS)
plot(NobsMOR,PRED.stat[2,],pch=16,ylim=range(c(NobsMOR,range(PRED.stat,na.rm=T))))
segments(x0=NobsMOR,x1=NobsMOR,y0=PRED.stat[1,],y1=PRED.stat[3,])
abline(a=0,b=1,lty=4,col=4,lwd=2)

#### fig 4

load("data/estim_effects.rda")

plot(x=0,y=0,type="n",ylim=c(0.5,12),xlim=c(-0.4,0.6))
for (i in 1:8){
  toto = Efaune[[i]][!is.na(Efaune[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,12),xlim=c(-0.05,0.1))
for (i in 1:8){
  toto = Evitesse[[i]][!is.na(Evitesse[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,12),xlim=c(-2,2))
for (i in 1:8){
  toto = Etrafic[[i]][!is.na(Etrafic[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,12),xlim=c(-0.2,0.5))
for (i in 1:8){
  toto = Edisthydro[[i]][!is.na(Edisthydro[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,12),xlim=c(-0.5,0.2))
for (i in 1:8){
  toto = Edistveg[[i]][!is.na(Edistveg[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

## alive guys ------------------------------------------------------------------

espece="Chevreuil europeen"

DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
datviv = DON[[1]]
datcoll = DON[[2]]
NobsVIV = DON[[3]]
NobsMOR = DON[[4]]

load(paste0("data/fitVIV_",espece,".rda"))
vivant=log(as.vector(predVIV[2,]))
v1<-ggplot()+gg(pxlVIV, aes(fill=vivant))+gg(datviv,size=0.5)
projgrid <- inla.mesh.projector(mesh, coordinates(pxlVIV))
Wmean <- inla.mesh.project(projgrid, ppVIV$summary.random$i$mean)
v2=ggplot()+gg(pxlVIV, aes(fill=Wmean))+gg(datviv,size=0.5)
v3=ggplot()+gg(pxlVIV, aes(fill=vivant-Wmean))+gg(datviv,size=0.5)

v1 # alive
v2 # only spatial field for alive guys
v3 # covariates (without spatial field)

## dead guys -------------------------------------------------------------------

load(paste0("data/pred_",espece,".rda"))

mort = apply(matrix(unlist(lapply(predMORT, function(x)
  as.vector(x))), ncol = 100), 1, mean)
m1 <- ggplot() + gg(pxlMOR, aes(fill = log(mort))) + gg(datcoll, size =
                                                          0.5)
m1 # prediction of dead guys
load(paste0("data/modINLA_",espece,"1.rda"))
projgrid <- inla.mesh.projector(meshROUTE, coordinates(pxlMOR))
Wmean <- inla.mesh.project(projgrid, modINLA$summary.random$i$mean)

v2=ggplot()+gg(pxlMOR, aes(fill=Wmean))+gg(datcoll,size=0.5)
v2 # spatial field of dead guys

library(patchwork)
m1 + v2


