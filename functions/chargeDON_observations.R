chargeDON = function(espece,Tmin,Tmax,CRSdon){
  load("data/donBON.rda")
  collision=roadkill_sp#read.csv(file="../data/collisions_mam_grand_non_carn.csv", h=T)
  collision=collision[(collision$annee>Tmin)&(collision$annee<Tmax),]#
  collision=collision[collision$spNEW==espece,]
  vivant=datVIV#read.csv(file="../data/observation_mam_grand_non_carn.csv", h=T)
  vivant=vivant[(vivant$annee>Tmin)&(vivant$annee<Tmax),]#
  vivant=vivant[vivant$nom_vern==espece,]
  
  covcol = data.frame(vitesse=f.vitesse(x=collision[,"X"]/1000,y=collision[,"Y"]/1000),
                      trafic=f.trafic(x=collision[,"X"]/1000,y=collision[,"Y"]/1000))#,collision[,c(1,4:6)])
  
  covviv = data.frame(sampling=f.sampling(x=vivant[,"x"]/1000,y=vivant[,"y"]/1000))#,vivant[,c(2,5)])

  covcol$trafic[covcol$trafic==0]=0.001

  datcoll <- SpatialPointsDataFrame(collision[,c("X", "Y")]/1000,covcol,proj4string=CRSdon)
  coordnames(datcoll) = c("x","y")
  datviv <- SpatialPointsDataFrame(vivant[,c("x", "y")]/1000,covviv,proj4string=CRSdon)
  coordnames(datviv) = c("x","y")
  
  #pour ajustement
  r = as(pxlVIV, "SpatialPolygonsDataFrame")
  rSF=st_as_sf(r)
  datvivSF = st_as_sf(datviv)
  test=st_intersects(rSF,datvivSF)
  NobsVIV=unlist(lapply(test,length))
  r = as(pxlMOR, "SpatialPolygonsDataFrame")
  rSF=st_as_sf(r)
  datcollSF = st_as_sf(datcoll)
  test=st_intersects(rSF,datcollSF)
  NobsMOR=unlist(lapply(test,length))
  
  return(list(datviv,datcoll,NobsVIV,NobsMOR))
}
