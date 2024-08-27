source("functions/inlabookfunctions.R")

####covariables, ACP et fonction pur aller chercher les covariables


foretR <- raster("data/raster_cov/foret_landcover_1km.tif")
extent(foretR) <- extent(c(xmin(foretR), xmax(foretR), ymin(foretR), ymax(foretR))/1000)
projection(foretR) <- gsub("units=m", "units=km", projection(foretR))
foretSG <- as(scale(foretR), 'SpatialGridDataFrame')
samplingR <- raster("data/raster_cov/sampling_axis1_1km.tif")
extent(samplingR) <- extent(c(xmin(samplingR), xmax(samplingR), ymin(samplingR), ymax(samplingR))/1000)
projection(samplingR) <- gsub("units=m", "units=km", projection(samplingR))
samplingSG <- as(-scale(exp(samplingR)), 'SpatialGridDataFrame')
cultureR <- raster("data/raster_cov/culture_landcover_1km.tif")
extent(cultureR) <- extent(c(xmin(cultureR), xmax(cultureR), ymin(cultureR), ymax(cultureR))/1000)
projection(cultureR) <- gsub("units=m", "units=km", projection(cultureR))
cultureSG <- as(scale(cultureR), 'SpatialGridDataFrame')
haieR <- raster("data/raster_cov/densite_haie_1km.tif")
extent(haieR) <- extent(c(xmin(haieR), xmax(haieR), ymin(haieR), ymax(haieR))/1000)
projection(haieR) <- gsub("units=m", "units=km", projection(haieR))
haieSG <- as(scale(haieR), 'SpatialGridDataFrame')
hydroR <- raster("data/raster_cov/densite_hydro_1km.tif")
extent(hydroR) <- extent(c(xmin(hydroR), xmax(hydroR), ymin(hydroR), ymax(hydroR))/1000)
projection(hydroR) <- gsub("units=m", "units=km", projection(hydroR))
hydroSG <- as(scale(hydroR), 'SpatialGridDataFrame')
urbainR <- raster("data/raster_cov/urbain_landcover_1km.tif")
extent(urbainR) <- extent(c(xmin(urbainR), xmax(urbainR), ymin(urbainR), ymax(urbainR))/1000)
projection(urbainR) <- gsub("units=m", "units=km", projection(urbainR))
urbainSG <- as(scale(urbainR), 'SpatialGridDataFrame')

vitesseR <- raster("data/raster_cov/vitesse_1km.tif")
extent(vitesseR) <- extent(c(xmin(vitesseR), xmax(vitesseR), ymin(vitesseR), ymax(vitesseR))/1000)
projection(vitesseR) <- gsub("units=m", "units=km", projection(vitesseR))
vitesseR[vitesseR==0] = NA
vitesseR = focal(vitesseR, w = matrix(1,11,11), fun = mean, 
                 pad = TRUE, na.rm = TRUE )
vitesseSG <- as(vitesseR, 'SpatialGridDataFrame')

traficR <- raster("data/raster_cov/trafic_1km.tif")
extent(traficR) <- extent(c(xmin(traficR), xmax(traficR), ymin(traficR), ymax(traficR))/1000)
projection(traficR) <- gsub("units=m", "units=km", projection(traficR))
traficR[traficR==0] = NA
traficR = focal(traficR, w = matrix(1,11,11), fun = mean, 
                pad = TRUE, na.rm = TRUE )
traficSG <- as(traficR, 'SpatialGridDataFrame')

hydrodistR <- raster("data/raster_cov/raster_hydro_dist_50m.tif")
extent(hydrodistR) <- extent(c(xmin(hydrodistR), xmax(hydrodistR), ymin(hydrodistR), ymax(hydrodistR))/1000)
projection(hydrodistR) <- gsub("units=m", "units=km", projection(hydrodistR))
hydrodistSG <- as(hydrodistR/1000, 'SpatialGridDataFrame')

haiedistR <- raster("data/raster_cov/raster_veg_dist_1km.tif")
extent(haiedistR) <- extent(c(xmin(haiedistR), xmax(haiedistR), ymin(haiedistR), ymax(haiedistR))/1000)
projection(haiedistR) <- gsub("units=m", "units=km", projection(haiedistR))
haiedistSG <- as(log(haiedistR+1), 'SpatialGridDataFrame')


covF=foretSG$foret_landcover_1km
covU=urbainSG$urbain_landcover_1km
covC=cultureSG$culture_landcover_1km
covHy=hydroSG$densite_hydro_1km
covH=haieSG$densite_haie_1km

datACP = data.frame(foret=covF,culture=covC,hydro=covHy,haie=covH,urbain=covU)#
rownames(datACP)=paste("L",1:nrow(datACP))
datACP = datACP[!is.na(apply(datACP,1,sum)),]

pc <- dudi.pca(datACP,scannf=F,nf=4,center = F, scale = F)


datACP=cbind(datACP,pc$li)
datRASTER = data.frame(foret=covF,urbain=covU,culture=covC,hydro=covHy,haie=covH)
rownames(datRASTER)=paste("L",1:nrow(datRASTER))
datNEW = matrix(NA,ncol=ncol(datACP),nrow=nrow(datRASTER))
linonNA = (1:nrow(datRASTER))[!is.na(apply(datRASTER,1,sum))]
datNEW[linonNA,] = as.matrix(datACP)
colnames(datNEW) = colnames(datACP)

Axis1SG=SpatialGridDataFrame(grid=foretSG@grid, data=data.frame(Axis1=datNEW[,"Axis1"]), proj4string = fm_sp_get_crs(foretSG))
Axis2SG=SpatialGridDataFrame(grid=foretSG@grid, data=data.frame(Axis2=datNEW[,"Axis2"]), proj4string = fm_sp_get_crs(foretSG))
Axis3SG=SpatialGridDataFrame(grid=foretSG@grid, data=data.frame(Axis3=datNEW[,"Axis3"]), proj4string = fm_sp_get_crs(foretSG))
Axis4SG=SpatialGridDataFrame(grid=foretSG@grid, data=data.frame(Axis4=datNEW[,"Axis4"]), proj4string = fm_sp_get_crs(foretSG))



f.Axis1 <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(Axis1SG))
  proj4string(spp) <- fm_sp_get_crs(Axis1SG)
  v <- over(spp, Axis1SG)
  v$Axis1[is.na(v$Axis1)]=0
  return(v$Axis1)
}

f.Axis2 <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(Axis2SG))
  proj4string(spp) <- fm_sp_get_crs(Axis2SG)
  v <- over(spp, Axis2SG)
  v$Axis2[is.na(v$Axis2)]=0
  return(v$Axis2)
}

f.Axis3 <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(Axis3SG))
  proj4string(spp) <- fm_sp_get_crs(Axis3SG)
  v <- over(spp, Axis3SG)
  v$Axis3[is.na(v$Axis3)]=0
  return(v$Axis3)
}

f.Axis4 <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(Axis4SG))
  proj4string(spp) <- fm_sp_get_crs(Axis4SG)
  v <- over(spp, Axis4SG)
  v$Axis4[is.na(v$Axis4)]=0
  return(v$Axis4)
}

f.foret <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(foretSG))
  proj4string(spp) <- fm_sp_get_crs(foretSG)
  v <- over(spp, foretSG)
  v$foret_landcover_1km[is.na(v$foret_landcover_1km)]=0
  return(v$foret_landcover_1km)
}

f.sampling <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(samplingSG))
  proj4string(spp) <- fm_sp_get_crs(samplingSG)
  v <- over(spp, samplingSG)
  v$layer[is.na(v$layer)]=0
  return(v$layer)
}

f.culture <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(cultureSG))
  proj4string(spp) <- fm_sp_get_crs(cultureSG)
  v <- over(spp, cultureSG)
  v$culture_landcover_1km[is.na(v$culture_landcover_1km)]=0
  return(v$culture_landcover_1km)
}

f.haie <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(haieSG))
  proj4string(spp) <- fm_sp_get_crs(haieSG)
  v <- over(spp, haieSG)
  v$densite_haie_1km[is.na(v$densite_haie_1km)]=0
  return(v$densite_haie_1km)
}

f.hydro <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(hydroSG))
  proj4string(spp) <- fm_sp_get_crs(hydroSG)
  v <- over(spp, hydroSG)
  if (any(is.na(v$densite_hydro_1km))) {
    v$densite_hydro_1km <- inlabru:::bru_fill_missing(hydroSG, spp, v$densite_hydro_1km)
  }
  return(v$densite_hydro_1km)
}

f.urbain <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(urbainSG))
  proj4string(spp) <- fm_sp_get_crs(urbainSG)
  v <- over(spp, urbainSG)
  v$urbain_landcover_1km[is.na(v$urbain_landcover_1km)]=0
  return(v$urbain_landcover_1km)
}


f.vitesse <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(vitesseSG))
  proj4string(spp) <- fm_sp_get_crs(vitesseSG)
  v <- over(spp, vitesseSG)
  v$layer[is.na(v$layer)]=0
  return(v$layer)
}

f.trafic <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(traficSG))
  proj4string(spp) <- fm_sp_get_crs(traficSG)
  v <- over(spp, traficSG)
  v$layer[is.na(v$layer)]=0
  return(v$layer)
}

f.disthydro <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(hydrodistSG))
  proj4string(spp) <- fm_sp_get_crs(hydrodistSG)
  v <- over(spp, hydrodistSG)
  v$raster_hydro_dist_50m[is.na(v$raster_hydro_dist_50m)]=0
  return(v$raster_hydro_dist_50m)
}

f.distveg <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(haiedistSG))
  proj4string(spp) <- fm_sp_get_crs(haiedistSG)
  v <- over(spp, haiedistSG)
  v$layer[is.na(v$layer)]=0
  return(v$layer)
}

###mesh2D
load("data/donBON.rda")
collision=roadkill_sp#read.csv(file="../data/collisions_mam_grand_non_carn.csv", h=T)
vivant=datVIV#read.csv(file="../data/observation_mam_grand_non_carn.csv", h=T)
domaine = st_read("data/domaine/emprise_maillage_1km_g2.shp")
domaine = st_cast(domaine,"POLYGON")
domaine = domaine[5,]
domaine=as_Spatial(domaine)
domaine=spTransform(domaine,"+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
domaine=st_as_sf(domaine)
coordDOM <- st_coordinates(domaine)
coordMESH=as.matrix(coordDOM[,1:2])
tutu=as.matrix(vivant[,c("x","y")]/1000)
coordMESH=rbind(coordMESH,tutu)
bndint = inla.nonconvex.hull(coordMESH, convex=-.03)
bndext = inla.nonconvex.hull(coordMESH, convex=-.2)
mesh = inla.mesh.2d(loc=rbind(as.matrix(collision[,c("X","Y")]),as.matrix(vivant[,c("x","y")]))/1000,
                    boundary = list(bndint,bndext),max.edge = c(3,15),cutoff = 3,
                    crs=proj4string(foretSG))

st_crs(domaine) <- NA
domaineSP <- as_Spatial(domaine)
proj4string(domaineSP)=proj4string(foretSG)


dmesh0 <- book.mesh.dual(mesh)
proj4string(dmesh0) <- proj4string(foretSG)
dmesh <- st_as_sf(dmesh0)
proj4string(domaineSP) <- proj4string(foretSG)
domaineSP1 <- st_as_sf(domaineSP)

w <- sapply(1:dim(dmesh)[1], function(i) {
  if (st_intersects(dmesh[i, ], domaineSP1,sparse = F))
    return(st_area(st_intersection(dmesh[i, ], domaineSP1)))
  else return(0)
})

#sum(w);area(domaineSP)
#table(w>0)
#par(mfrow=c(1,1))
#plot(dmesh,col=(w>0)*4)

###mesh routes
routes = st_read("data/routesMAELISS/routes_diro_fill.shp")
routes=as_Spatial(routes)
routes=spTransform(routes,"+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
routes=st_as_sf(routes)
routes <- st_crop(routes,c(xmin=extent(domaine)[1],ymin=extent(domaine)[3],xmax=extent(domaine)[2],ymax=extent(domaine)[4]))

routeBUFtmp=st_buffer(routes,dist=2)
routeBUFtmp=st_union(routeBUFtmp)

routesSL <- as(routeBUFtmp,"Spatial")
proj4string(routesSL)=proj4string(foretSG)

routeBUF0=st_buffer(routeBUFtmp,dist=-1)
routeBUF0 <- as(routeBUF0,"Spatial")
routeBUF <- st_as_sf(routeBUF0)
proj4string(routeBUF)=mesh$crs
routeBUF <- st_transform(routeBUF, mesh$crs)

meshROUTE = inla.mesh.2d(boundary = bndint,max.edge = 3,cutoff = 0.5,crs=mesh$crs)
route.tri = inla.over_sp_mesh(routeBUF, y = meshROUTE, 
                              type = "centroid", ignore.CRS = TRUE)
num.tri <- length(meshROUTE$graph$tv[, 1])
barrier.tri <- setdiff(1:num.tri, route.tri)
poly.barrier <- inla.barrier.polygon(meshROUTE, 
                                     barrier.triangles = barrier.tri)

dmeshROUTE0 <- book.mesh.dual(meshROUTE)
dmeshROUTE <- st_as_sf(dmeshROUTE0)
st_crs(dmeshROUTE) <- mesh$crs


wROUTE <- sapply(1:length(dmeshROUTE), function(i) {
  if (st_intersects(dmeshROUTE[i, ], routeBUF,sparse = F))
    return(st_area(st_intersection(dmeshROUTE[i, ], routeBUF)))
  else return(0)
  })
#plot(dmeshROUTE,col=(wROUTE>0)*4)

## pour ajustement
toto = extent(domaineSP)
xrange = toto[2]-toto[1]
yrange = toto[4]-toto[3]
resopred=80
pxlVIV <- fm_pixels(mesh, mask=domaineSP, dims = c(resopred,  round(resopred*yrange/xrange)))
projgridVIV <- inla.mesh.projector(mesh, st_coordinates(pxlVIV))
resopred=300
pxlMOR <- fm_pixels(meshROUTE, mask=routeBUF, dims = c(resopred,  round(resopred*yrange/xrange)))
projgridMOR <- inla.mesh.projector(meshROUTE, st_coordinates(pxlMOR))

save.image("OUT_chargeDON_covariables_mesh_modifJYB31072024.rda")





