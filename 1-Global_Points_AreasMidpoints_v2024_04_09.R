# Load the required libraries
library(terra)         # CRAN v1.7-29
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0



# Check the zip file, get species names -----------------------------------

#process points #####
points <- unzip("temp/points_v1.zip",list = T)

ptifs <- grep('\\.tif',points$Name)
ptifs <- points[ptifs,]

#select only the species names

ptifs_sp <- gsub('.*/X__',"",ptifs$Name)
ptifs_sp <- gsub('__points',"",ptifs_sp)
ptifs_sp <- gsub('\\.tif',"",ptifs_sp)
ptifs_sp <- data.frame(names=ptifs_sp)
ptifs_sp <- ptifs_sp$names


# Prepare objects to be used below ----------------------------------------

w <- vect("wrld_moll.shp")


# get area and unzip ------------------------------------------------------

#points
plist <- ptifs$Name

resu_p <- tibble(midX=NA ,midY=NA,midX84=NA, midY84=NA,species=NA, areaM=NA,area84=NA,region=NA,band=NA,continent=NA,con2=NA,con3=NA,con4=NA,con5=NA)

#for(j in 1:length(ptifs$Name)){
  for(j in 1:150){
  print(j)
  
ras <- rast(unzip("temp/points_v1.zip",files = plist[j] ))
resu_p[j,"species"] <- ptifs_sp[j]
  
# Calculate midpoints -----------------------------------------------------
p <- terra::as.polygons(ras)  
p84 <- project(p,"epsg:4326")
  
c <- terra::centroids(p,inside=FALSE)
cc <- terra::geom(c)[, 3:4, drop = FALSE]
  
resu_p[j,"midX"] <- cc[1]
resu_p[j,"midY"] <- cc[2]
  
c84 <- terra::centroids(p84,inside=FALSE)
cc84 <- terra::geom(c84)[, 3:4, drop = FALSE]

resu_p[j,"midX84"] <- cc84[1]
resu_p[j,"midY84"] <- cc84[2]

# Classify species according to their area position -----------------------
  
ymin <- as.vector(terra::ext(p)[3])
ymax <- as.vector(terra::ext(p)[4])
  
xmin <- as.vector(terra::ext(p)[1])
xmax <- as.vector(terra::ext(p)[2])
  

if(F){ # to calculate using lat/long
region <- NA
  if(ymax < -23.4 | ymin > 23.4) {region <- "Temperate"}
  if(ymin > -23.4 & ymax < 23.4) {region <- "Tropical"}
  if(is.na(region))              {region <- "Temperate/Tropical"}
  
  resu[j,"region"] <- region
  
  band <- NA
  if(xmin > -169 & xmax  < -26)  {band <- "Americas"}
  if(xmin > -26  & xmax  <  52)  {band <- "Africa/Europe"}
  if(xmin > 52   & xmax  < 179) {band <- "Oceania/Asia"}
  
  resu[j,"band"] <- band
}

if(T){ # to calculate using lat/long Mollweide
  region <- NA
  if(ymax < -2862317 | ymin > 2862317) {region <- "Temperate"}
  if(ymin > -2862317 & ymax < 2862317) {region <- "Tropical"}
  if(is.na(region))              {region <- "Temperate/Tropical"}
  
  resu_p[j,"region"] <- region
  
  band <- NA
  if(xmin > -15280020.4 & xmax  < -364445.7)  {band <- "Americas"}
  if(xmin > -1730360    & xmax  <  6000000 )  {band <- "Africa/Europe"}
  if(xmin >  6000000     & xmax  < 17601618 ) {band <- "Oceania/Asia"}
  
  resu_p[j,"band"] <- band
}

  
# get continents ----------------------------------------------------------
  
  f <- crop(w,p)
  co<- f$region_un
  
  if(length(unique(co)) == 0){resu_p[j,"continent"] <- NA } 
  if(length(unique(co)) == 1){resu_p[j,"continent"] <- co[1] } 
  if(length(unique(co)) >1){
    for(u in 1:length(unique(co))){
      nn <- which(colnames(resu_p)=="continent")
      resu_p[j,nn-1+u] <- unique(co)[u] }}
  
  
# Calculate range areas ---------------------------------------------------
  
area <- expanse(p,unit="km")
resu_p[j,"areaM"] <- area
  
area84 <- expanse(p84,unit="km")
resu_p[j,"area84"] <- area84

}

write_rds(resu_p,"products/areas_midpoints_global_moll/area_midpoints_POINTS_regionsGLOBAL.rds")


