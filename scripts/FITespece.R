#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(sf)
library(inlabru)
library(INLA)
library(ggplot2)
library(maptools)
library(igraph)
library(dplyr)
library(raster)
library(tidygraph)
library(tidyverse)
library(ade4)
library(fields)
library(rgeos)
library(pROC)

setwd("BON")

#la fonction pixels d'inlabru ne marche pas toujours...
#source("pixelBIS.R")

#covariables et mesh
source("chargeDON_covariables_mesh.R")

#observations
source("chargeDON_observations.R")
listeSP=intersect(unique(roadkill_sp$spNEW),unique(datVIV$nom_vern))
for (espece in listeSP[6]){
DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
datviv = DON[[1]]
datcoll = DON[[2]]
NobsVIV = DON[[3]]
NobsMOR = DON[[4]]

png(file=paste0("don_",espece,".png"),width=300,height=200,res=300,units="mm")
par(mar=rep(1,4),mfrow=c(1,1))
plot(mesh, main = "", asp = 1)
plot(routesSL,add=T,col=2,border=2,lwd=2)
plot(domaineSP,add=T,border=1,lwd=2)
plot(datviv,add=T,col=3,pch=16,cex=0.5)
plot(datcoll,add=T,col=4,pch=16,cex=0.5)
dev.off()

################
## Faune vivante
################


CAxis1 = f.Axis1(x=c(mesh$loc[,1],coordinates(datviv)[,1]),
                 y=c(mesh$loc[,2],coordinates(datviv)[,2]))
CAxis2 = f.Axis2(x=c(mesh$loc[,1],coordinates(datviv)[,1]),
                 y=c(mesh$loc[,2],coordinates(datviv)[,2]))
CAxis3 = f.Axis3(x=c(mesh$loc[,1],coordinates(datviv)[,1]),
                 y=c(mesh$loc[,2],coordinates(datviv)[,2]))
CAxis4 = f.Axis4(x=c(mesh$loc[,1],coordinates(datviv)[,1]),
                 y=c(mesh$loc[,2],coordinates(datviv)[,2]))
Csampling = f.sampling(x=c(mesh$loc[,1],coordinates(datviv)[,1]),
                       y=c(mesh$loc[,2],coordinates(datviv)[,2]))

### estimation du range
matern <- inla.spde2.pcmatern(mesh,
                              alpha = 2, # fractional operator which is related 
                              #to the smoothness of the Gaussian field
                              prior.sigma = c(1.6, 0.1), # P(sigma > 1) = 0.5
                              prior.range = c(55, 0.9)) # P(range < 100) = 0.9

nv = mesh$n
n = nrow(datviv)
y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))
imat <- Diagonal(nv, rep(1, nv))
lmat <- inla.spde.make.A(mesh, coordinates(datviv))

A.pp <- rbind(imat, lmat)
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp),
  A = list(1, A.pp), 
  effects = list(list(b0 = 1, Axis1=CAxis1, Axis2=CAxis2, Axis3=CAxis3, Axis4=CAxis4, sampling=Csampling),
                 list(i = 1:nv)),
  tag = 'pp')

ppVIV <- inla(y ~ 0 + b0 + f(i, model = matern),#Axis1 + Axis2+ Axis3 + Axis4 + sampling,#
              family = 'poisson', data = inla.stack.data(stk.pp), 
              control.predictor = list(A = inla.stack.A(stk.pp)),
              E = inla.stack.data(stk.pp)$e,
              control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE, return.marginals.predictor=TRUE), 
              control.inla=list(strategy="simplified.laplace",int.strategy="eb"))

sumsum=summary(ppVIV)
rangeFIX=sumsum$hyperpar$mean[1]
sigmaFIX=sumsum$hyperpar$mean[2]

## estimation des effets avec le range fixé
matern <- inla.spde2.pcmatern(mesh,
                              alpha = 2, # fractional operator which is related 
                              #to the smoothness of the Gaussian field
                              prior.sigma = c(sigmaFIX, NA), # P(sigma > 1) = 0.5
                              prior.range = c(rangeFIX, NA)) # P(range < 100) = 0.9


ppVIV <- inla(y ~ 0 + b0  + Axis1 + Axis2+ Axis3 + Axis4 + sampling + f(i, model = matern),#,#
              family = 'poisson', data = inla.stack.data(stk.pp), 
              control.predictor = list(A = inla.stack.A(stk.pp)),
              E = inla.stack.data(stk.pp)$e,
              control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE, return.marginals.predictor=TRUE), 
              control.inla=list(strategy="simplified.laplace",int.strategy="eb"))

fitvalmin = inla.mesh.project(projgridVIV, ppVIV$summary.fitted.values$'0.025quant'[1:mesh$n])
fitvalmed = inla.mesh.project(projgridVIV, ppVIV$summary.fitted.values$'0.5quant'[1:mesh$n])
fitvalmax = inla.mesh.project(projgridVIV, ppVIV$summary.fitted.values$'0.975quant'[1:mesh$n])
r = as(pxlVIV, "SpatialPolygonsDataFrame")
rSF=st_as_sf(r)
predVIV=rbind(fitvalmin,fitvalmed,fitvalmax)*st_area(rSF)
save(ppVIV,predVIV,NobsVIV,file=paste0("fitVIV_",espece,".rda"))


## prédiction d'inla en raster pour ensuite l'injecter dans le modèle route, estimation du modèle route

Nrep=100
pr.int.tot <- inla.posterior.sample(Nrep,ppVIV)

RES=list()
predMORT=list()
predAUC=list()
for (i in 1:Nrep){
  predTMP = inla.mesh.project(projgridVIV, pr.int.tot[[i]]$latent[1:mesh$n])
  predTMP_R = pxlVIV
  predTMP_R@data=as.data.frame(predTMP)
  FvivSG <- as(predTMP_R, 'SpatialGridDataFrame')
  
  f.Fviv <- function(x, y) {
    spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(FvivSG))
    proj4string(spp) <- fm_sp_get_crs(FvivSG)
    v <- over(spp, FvivSG)
    v$predTMP[is.na(v$predTMP)]=0
    return(v$predTMP)
  }
  
  Nbebete = f.Fviv(x=c(meshROUTE$loc[,1],coordinates(datcoll)[,1]),
                   y=c(meshROUTE$loc[,2],coordinates(datcoll)[,2]))
  vitesse = f.vitesse(x=c(meshROUTE$loc[,1],coordinates(datcoll)[,1]),
                      y=c(meshROUTE$loc[,2],coordinates(datcoll)[,2]))
  trafic = f.trafic(x=c(meshROUTE$loc[,1],coordinates(datcoll)[,1]),
                    y=c(meshROUTE$loc[,2],coordinates(datcoll)[,2]))
  
  maternROUTE <- inla.barrier.pcmatern(meshROUTE,barrier.triangles = barrier.tri,
                                       prior.sigma = c(1, 0.5),#NA
                                       prior.range = c(15, 0.5))
  
  
  nv = meshROUTE$n
  n = nrow(datcoll)
  y.pp <- rep(0:1, c(nv, n))
  e.pp <- c(wROUTE, rep(0, n))
  imat <- Diagonal(nv, rep(1, nv))
  lmat <- inla.spde.make.A(meshROUTE, coordinates(datcoll))
  
  A.pp <- rbind(imat, lmat)
  stk.pp <- inla.stack(
    data = list(y = y.pp, e = e.pp),
    A = list(1, A.pp), 
    effects = list(list(b0 = 1, Nbebete = Nbebete,vitesse=vitesse, trafic=log(trafic+1)),
                   list(i = 1:nv)),
    tag = 'pp')
  
  ppMOR <- inla(y ~ 0 + b0 + Nbebete + vitesse + trafic + f(i, model = maternROUTE),
                   family = 'poisson', data = inla.stack.data(stk.pp), 
                   control.predictor = list(A = inla.stack.A(stk.pp)),
                   E = inla.stack.data(stk.pp)$e,
                   control.compute=list(dic=TRUE, cpo=TRUE, config = TRUE, return.marginals.predictor=TRUE), 
                   control.inla=list(strategy="simplified.laplace",
                                     int.strategy="eb"))
  
  RES[[i]]=summary(ppMOR)
  
  fitval = inla.mesh.project(projgridMOR, ppMOR$summary.fitted.values$'0.5quant'[1:meshROUTE$n])
  r = as(pxlMOR, "SpatialPolygonsDataFrame")
  rSF=st_as_sf(r)
  predMORT[[i]]=fitval*st_area(rSF)
  predAUC[[i]] = 1-exp(-as.numeric(fitval*st_area(rSF)))
}
save(RES,predMORT,predAUC,NobsMOR, file=paste0("RES_",espece,".rda"))
}


##quelques graphes
source("chargeDON_observations.R")
load("res/OUT_chargeDON_covariables_mesh.rda")
listeSP=intersect(unique(roadkill_sp$spNEW),unique(datVIV$nom_vern))
listeSP
espece="Sanglier"

DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
datviv = DON[[1]]
datcoll = DON[[2]]
NobsVIV = DON[[3]]
NobsMOR = DON[[4]]

load(paste0("res/fitVIV_",espece,".rda"))
load(paste0("res/RES_",espece,".rda"))

Nrep=length(RES)
NobsMOR_bis = as.numeric(NobsMOR > 0) # les cellules où il y a présence
AUC=NULL
TSS=NULL
for (i in 1:Nrep){
  testroc=roc(NobsMOR_bis, predAUC[[i]])
  AUC=c(AUC,testroc$auc)
  TSS=c(TSS,max(testroc$sensitivities+testroc$specificities-1))
}
par(mfrow=c(1,2))
hist(AUC)
hist(TSS)

#testroc=roc(NobsMOR_bis, predAUC[[1]])
#testroc # AUC
#plot.roc(NobsMOR_bis, predAUC[[1]], col="red", print.thres="best")
#AUC=testroc$auc
#TSS=testroc$sensitivities+testroc$specificities-1
#plot(testroc$thresholds,TSS,type="l")



predMM = matrix(unlist(predMORT),ncol=length(predMORT[[1]]),byrow=T)
predMM_stat = apply(predMM,2,quantile,probs=c(0.025,0.5,0.975))
par(mfrow=c(1,1))
plot(NobsMOR,predMM_stat[2,],pch=16,ylim=range(predMM_stat,na.rm=T))
segments(x0=NobsMOR,x1=NobsMOR,y0=predMM_stat[1,],y1=predMM_stat[3,])
abline(a=0,b=1,lty=4,col=4,lwd=2)


Efaune=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["Nbebete",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
plot(1:100,Efaune[,2],pch=16,ylim=range(Efaune,na.rm=T))
segments(x0=1:100,x1=1:100,y0=Efaune[,1],y1=Efaune[,3])
abline(h=mean(Efaune[,2],na.rm=T),lty=4,col=4,lwd=2)

Evitesse=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["vitesse",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
plot(1:100,Evitesse[,2],pch=16,ylim=range(Evitesse,na.rm=T))
segments(x0=1:100,x1=1:100,y0=Evitesse[,1],y1=Evitesse[,3])
abline(h=mean(Evitesse[,2],na.rm=T),lty=4,col=4,lwd=2)

Etrafic=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["trafic",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
plot(1:100,Etrafic[,2],pch=16,ylim=range(Etrafic,na.rm=T))
segments(x0=1:100,x1=1:100,y0=Etrafic[,1],y1=Etrafic[,3])
abline(h=mean(Etrafic[,2],na.rm=T),lty=4,col=4,lwd=2)



projgrid <- inla.mesh.projector(mesh, coordinates(pxlVIV))
Wmean <- inla.mesh.project(projgrid, ppVIV$summary.random$i$mean)
ggplot()+gg(pxlVIV, aes(fill=Wmean))+gg(datviv,size=0.5)

ggplot()+gg(pxlVIV, aes(fill=log(as.vector(predVIV[2,]))))+gg(datviv,size=0.5)

ggplot()+gg(pxlMOR, aes(fill=log(as.vector(predMORT[[1]]))))+gg(datcoll,size=0.5)


#marfitted <- ppMOR$marginals.fitted.values[1:meshROUTE$n]
#threshold <- 0.4
#exceed.prob <- lapply(X= marfitted, FUN = function(x) inla.pmarginal(marginal = x, threshold))
#exceed.prob <- 1 - unlist(exceed.prob)
#exceed.prob <- inla.mesh.project(projgrid, exceed.prob)
#ggplot()+gg(pxl, aes(fill=exceed.prob))+gg(datcoll,size=0.5)
#ggplot()+gg(pxlMOR, aes(fill=fitval))+gg(datcoll,size=0.5)
#plot(NobsMOR,fitval*st_area(rSF))
#abline(a=0,b=1,col=4,lwd=2,lty=2)





###fig 1 = mesh + données

###fig 3

##quelques graphes
source("chargeDON_observations.R")
load("res/OUT_chargeDON_covariables_mesh.rda")
listeSP=intersect(unique(roadkill_sp$spNEW),unique(datVIV$nom_vern)) [-6]
listeSP

AUCtot=NULL
TSStot=NULL
NobsMORtot=NULL
predSTATtot=list()
il=0
for (espece in listeSP){
il=il+1
DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
datviv = DON[[1]]
datcoll = DON[[2]]
NobsVIV = DON[[3]]
NobsMOR = DON[[4]]

load(paste0("res/fitVIV_",espece,".rda"))
load(paste0("res/RES_",espece,".rda"))

Nrep=length(RES)
NobsMOR_bis = as.numeric(NobsMOR > 0) # les cellules où il y a présence
AUC=NULL
TSS=NULL
for (i in 1:Nrep){
  testroc=roc(NobsMOR_bis, predAUC[[i]])
  AUC=c(AUC,testroc$auc)
  TSS=c(TSS,max(testroc$sensitivities+testroc$specificities-1))
}

predMM = matrix(unlist(predMORT),ncol=length(predMORT[[1]]),byrow=T)
predMM_stat = apply(predMM,2,quantile,probs=c(0.025,0.5,0.975))

AUCtot=cbind(AUCtot,AUC)
TSStot=cbind(TSStot,TSS)
NobsMORtot=cbind(NobsMORtot,NobsMOR)
predSTATtot[[il]]=predMM_stat

}

png(file="fit.png",units="mm",width=150,height=300,res=200)
par(mfrow=c(5,2))
for (i in 1:10){
  plot(NobsMORtot[,i],predSTATtot[[i]][2,],pch=16,ylim=range(c(NobsMORtot[,i],range(predSTATtot[[i]],na.rm=T))),
       main=paste(listeSP[i],"AUC=",round(mean(AUCtot[,i]),2),"+/-",round(sd(AUCtot[,i]),2),"TSS=",round(mean(TSStot[,i]),2),"+/-",round(sd(TSStot[,i]),2)))
  segments(x0=NobsMORtot[,i],x1=NobsMORtot[,i],y0=predSTATtot[[i]][1,],y1=predSTATtot[[i]][3,])
  abline(a=0,b=1,lty=4,col=4,lwd=2)
}
dev.off()


#### fig 4


Efaune=list()
Evitesse=list()
Etrafic=list()
il=0
for (espece in listeSP){
  il=il+1
  
  DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
  datviv = DON[[1]]
  datcoll = DON[[2]]
  NobsVIV = DON[[3]]
  NobsMOR = DON[[4]]
  
  load(paste0("res/fitVIV_",espece,".rda"))
  load(paste0("res/RES_",espece,".rda"))
  
  Efaunetmp=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["Nbebete",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
  Evitessetmp=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["vitesse",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
  Etrafictmp=matrix(unlist(lapply(RES,function(x){tutu=x$fixed["trafic",3:5];if(is.null(tutu)){tutu=rep(NA,3)};tutu})),ncol=3,byrow=T)
  
  Efaune[[il]]=Efaunetmp
  Evitesse[[il]]=Evitessetmp
  Etrafic[[il]]=Etrafictmp
}

plot(x=0,y=0,type="n",ylim=c(0.5,14),xlim=c(-0.4,0.6))
for (i in 1:10){
toto = Efaune[[i]][!is.na(Efaune[[i]][,2]),]
polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,14),xlim=c(-0.05,0.1))
for (i in 1:10){
  toto = Evitesse[[i]][!is.na(Evitesse[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)

plot(x=0,y=0,type="n",ylim=c(0.5,14),xlim=c(-2,2))
for (i in 1:10){
  toto = Etrafic[[i]][!is.na(Etrafic[[i]][,2]),]
  polygon(x=c(toto[,1],rev(toto[,3])),y=c(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),rev(seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)))),col="goldenrod",border="goldenrod")
  points(toto[,2],y=seq(i+0.3*(i-1),i+1+0.3*(i-1),length=nrow(toto)),pch=16,cex=0.5)
}
abline(v=0,col="steelblue",lty="dashed",lwd=2)


###fig6
espece="Chevreuil europeen"

DON = chargeDON(espece,2015,2020,CRSdon=CRS(proj4string(foretSG)))
datviv = DON[[1]]
datcoll = DON[[2]]
NobsVIV = DON[[3]]
NobsMOR = DON[[4]]

load(paste0("res/fitVIV_",espece,".rda"))
load(paste0("res/RES_",espece,".rda"))

library(patchwork)
vivant=log(as.vector(predVIV[2,]))
v1<-ggplot()+gg(pxlVIV, aes(fill=vivant))+gg(datviv,size=0.5)
projgrid <- inla.mesh.projector(mesh, coordinates(pxlVIV))
Wmean <- inla.mesh.project(projgrid, ppVIV$summary.random$i$mean)
v2=ggplot()+gg(pxlVIV, aes(fill=Wmean))+gg(datviv,size=0.5)
v3=ggplot()+gg(pxlVIV, aes(fill=vivant-Wmean))+gg(datviv,size=0.5)

mort=apply(matrix(unlist(lapply(predMORT,function(x) as.vector(x))),ncol=100),1,mean)
m1<-ggplot()+gg(pxlMOR, aes(fill=mort))+gg(datcoll,size=0.5)










