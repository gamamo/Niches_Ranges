# Load the required libraries
library(terra)         # CRAN v1.7-29
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0
library(future)
library(doFuture)
library(tictoc)


# Check the zip file, get species names -----------------------------------

#process range bagging #####
rb <- unzip("temp/rangebag_v1.zip", list=T)

rbt <- grep('\\.tif'  ,rb$Name)
rbtifs <- rb[rbt, ]
x016   <- grep("/X0.165_",rbtifs$Name)
rbtifs <- rbtifs[x016, ]
occ <- grep("_occEco.tif",rbtifs$Name)
rbtifs <- rbtifs[occ, ]

# get only species names
rbtifs_sp <- gsub("/X0.1.*","",rbtifs$Name)  
rbtifs_sp <- gsub(".*ps/","",rbtifs_sp)
rbtifs_sp <- data.frame(names=rbtifs_sp)
rbtifs_sp <- rbtifs_sp$names


# Prepare objects to be used below ----------------------------------------

w <- vect("wrld_moll.shp")

# get area and unzip ------------------------------------------------------

#points
plist <- rbtifs$Name

#split the data
num_parts <- 20

# Use cut to create a factor with three levels, representing the parts
groups <- cut(seq_along(plist), breaks = num_parts, labels = FALSE)
groups_sp <- cut(seq_along(rbtifs_sp), breaks = num_parts, labels = FALSE)

# Use split to divide the vector into three parts based on the groups factor
split_vector <- split(plist, groups)
split_vector_sp <- split(rbtifs_sp, groups)

plan(multisession)

for(n in 1:length(split_vector)){
  print(n)

#select which will be group
rblist <- split_vector[[n]]
rbtifs_sp <- split_vector_sp[[n]]

resu_p <- tibble(midX=NA ,midY=NA,midX84=NA, midY84=NA,species=NA, areaM=NA,area84=NA,region=NA,band=NA,continent=NA,con2=NA,con3=NA,con4=NA,con5=NA)


tic()
ff <- foreach(j=1:length(rblist), .combine = rbind) %dofuture% {
 # ff <- foreach(j=1:2, .combine = rbind) %dofuture% {
  
  #print(j)
  
  ras <- rast(unzip("temp/rangebag_v1.zip",files = plist[j] ))
  #resu_p[j,"species"] <- rbtifs_sp[j]
  species <- rbtifs_sp[j]
  
  # Calculate midpoints -----------------------------------------------------
  p <- terra::as.polygons(ras)  
  p84 <- project(p,"epsg:4326")
  
  c <- terra::centroids(p,inside=FALSE)
  cc <- terra::geom(c)[, 3:4, drop = FALSE]
  
  #resu_p[j,"midX"] <- cc[1]
  #resu_p[j,"midY"] <- cc[2]
  
  midX <- cc[1]
  midY <- cc[2]
  
  c84 <- terra::centroids(p84,inside=FALSE)
  cc84 <- terra::geom(c84)[, 3:4, drop = FALSE]
  
  #resu_p[j,"midX84"] <- cc84[1]
  #resu_p[j,"midY84"] <- cc84[2]
  
  midX84 <- cc84[1]
  midY84 <- cc84[2]
  
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
    
    #resu_p[j,"region"] <- region
    
    band <- NA
    if(xmin > -15280020.4 & xmax  < -364445.7)  {band <- "Americas"}
    if(xmin > -1730360    & xmax  <  6000000 )  {band <- "Africa/Europe"}
    if(xmin >  6000000     & xmax  < 17601618 ) {band <- "Oceania/Asia"}
    
    #resu_p[j,"band"] <- band
  }
  
  
  # get continents ----------------------------------------------------------
  if(F){
  f <- crop(w,p)
  co<- f$region_un
  
  if(length(unique(co)) == 0){resu_p[j,"continent"] <- NA } 
  if(length(unique(co)) == 1){resu_p[j,"continent"] <- co[1] } 
  if(length(unique(co)) >1){
    for(u in 1:length(unique(co))){
      nn <- which(colnames(resu_p)=="continent")
      resu_p[j,nn-1+u] <- unique(co)[u] }}
  }
  
  # Calculate range areas ---------------------------------------------------
  
  area <- expanse(p,unit="km")
  #resu_p[j,"areaM"] <- area
  
  area84 <- expanse(p84,unit="km")
  #resu_p[j,"area84"] <- area84
  
  c(species, midX, midY,midX84,midY84,area,area84, region, band)

}
toc()

ff2 <- as_tibble(ff)
ff3 <- ff2 |> rename(species=V1,midX=V2, midY=V3, midX84=V4,midY84=V5, area=V6,area84=V7,region=V8,band=V9) |>  
  mutate_at(2:6, as.numeric)

write_rds(ff3,
          paste0("products/areas_midpoints_global_moll/corrigido/area_midpoints_RANGEBAGGING_regionsGLOBAL",
                "_",n,"_",
                ".rds"))

}

ff3

