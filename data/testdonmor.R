roadkill = read.csv("C:/Users/jpapaix.DESKTOP-VU8SI74/Desktop/data/data_diro/roadkill_diro.csv")
head(roadkill)
unique(roadkill$espece)
range(roadkill$annee)

roadkill$espece[roadkill$espece=="Autre petit must\xe9lid\xe9"]="Autre petit mustelide"
roadkill$espece[roadkill$espece=="Li\xe8vre"]="Lievre"
roadkill$espece[roadkill$espece=="Rat musqu\xe9"]="Rat musque"
roadkill$espece[roadkill$espece=="\xc9cureuil"]="ecureuil"
roadkill$espece[roadkill$espece=="H\xe9risson"]="Herisson"
roadkill$espece[roadkill$espece=="Vison d'Am\xe9rique"]="Vison d'Amerique"

roadkill$obs[roadkill$obs==""]=NA


Lsp=read.csv("listesp.csv",sep=";",header=F)
Lsp=Lsp[,1:2]
colnames(Lsp)=c("Nsp","Gsp")

roadkill$espece[roadkill$espece=="Sangliers"]="Sanglier"
roadkill$espece[roadkill$espece=="Chevreuil"]="Chevreuil europeen"
roadkill$espece[roadkill$espece=="Vison d'Amerique"]="Vison damerique"
roadkill$espece[roadkill$espece=="Renard"]="Renard roux"
roadkill$espece[roadkill$espece=="ecureuil"]="Ecureuil"
roadkill$espece[roadkill$espece=="Putois"]="Putois deurope"
roadkill$espece[roadkill$espece=="Martre"]="Martre des pins"
roadkill$espece[roadkill$espece=="Blaireau"]="Blaireau europeen"
roadkill$espece[roadkill$espece=="Lievre"]="Lievre deurope"
roadkill$espece[roadkill$espece=="Amphibiens"]="Amphibien"
roadkill$espece[roadkill$espece=="Belette"]="Belette deurope"
roadkill$espece[roadkill$espece=="Cerf/Biche"]="Cerf elaphe"

roadkill$espece

spNEW = roadkill$espece
spNEW[spNEW=="Autre oiseau" & !is.na(spNEW)] = NA
spNEW[spNEW=="Rapaces diurnes" & !is.na(spNEW)] = NA
spNEW[spNEW=="Autre petite faune" & !is.na(spNEW)] = NA
spNEW[spNEW=="Reptiles" & !is.na(spNEW)] = NA
spNEW[spNEW=="Rapaces nocturnes" & !is.na(spNEW)] = NA
spNEW[spNEW=="Autre petit mustelide" & !is.na(spNEW)] = NA
spNEW[spNEW=="Autre rongeur" & !is.na(spNEW)] = NA
spNEW[spNEW=="Amphibien" & !is.na(spNEW)] = NA
spNEW[spNEW=="Chauve-souris" & !is.na(spNEW)] = NA
unique(spNEW)

grNEW = roadkill$espece
grNEW[grNEW=="Sanglier" | grNEW=="Cerf elaphe" | grNEW=="Chevreuil europeen"] = "Ongules"
grNEW[grNEW=="Hermine" | grNEW=="Fouine" | grNEW=="Martre des pins" | grNEW=="Belette deurope" | grNEW=="Putois deurope"] = "Mustelides"
grNEW[grNEW=="Blaireau europeen" | grNEW=="Renard roux" | grNEW=="Rat"] = "Mammiferes charognards"
grNEW[grNEW=="Lievre deurope" | grNEW=="Lapin" | grNEW=="Ecureuil" | grNEW=="Herisson"] = "Autres mammiferes"
grNEW[grNEW=="Rat musque" | grNEW=="Ragondin" | grNEW=="Loutre" | grNEW=="Vison damerique" | grNEW=="Castor"] = "Mammiferes aquatiques"
grNEW[grNEW=="Autre petit mustelide"] = "Mustelides"
grNEW[grNEW=="Autre rongeur"] = "Autres mammiferes"
unique(grNEW)

roadkill$spNEW=spNEW
roadkill$grNEW=grNEW
head(roadkill)

roadkill_sp=roadkill[!is.na(roadkill$spNEW),]

unique(roadkill_sp$spNEW)


###

datVIV=read.csv("C:/Users/jpapaix.DESKTOP-VU8SI74/Desktop/data/data_visionature/DFmam_format_final.csv")
unique(datVIV$nom_vern)

datVIV$nom_vern[datVIV$nom_vern=="Belette d'Europe"]="Belette deurope"
datVIV$nom_vern[datVIV$nom_vern=="Chevreuil européen"]="Chevreuil europeen"
datVIV$nom_vern[datVIV$nom_vern=="Lièvre d'Europe"]="Lievre deurope"
datVIV$nom_vern[datVIV$nom_vern=="Putois d'Europe"]="Putois deurope"
datVIV$nom_vern[datVIV$nom_vern=="Blaireau européen"]="Blaireau europeen"
datVIV$nom_vern[datVIV$nom_vern=="Cerf élaphe"]="Cerf elaphe"

save(datVIV,roadkill_sp,file="donBON.rda")









