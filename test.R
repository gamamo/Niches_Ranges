library(BIEN)
dir.create(path = "temp/BIEN2_ranges",recursive = TRUE)
BIEN_ranges_species_bulk(species = NULL,
                         directory = "temp/BIEN2_ranges/",
                         batch_size = 100,
                         use_parallel = FALSE)


# try the letsR -----------------------------------------------------------
library(letsR)
library(raster)
library(terra)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(phyloregion)

data(Phyllomedusa)
class(Phyllomedusa)
Phyllomedusa[[1]]

PAM <- lets.presab(Phyllomedusa, xmn = -93, xmx = -29,
                   ymn = -57, ymx = 15)
head(PAM)

fdir <- system.file("NGAplants", package="phyloregion")
s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
plot(s)
sp <- random_species(100, species=5, shp=s)
plot(sp)
class(sp)
plot(s,add=T)
pol <- polys2comm(dat = sp, species = "species")
head(pol[[1]])



am <-  st_read("temp/BIEN2_ranges/2/Abolboda_americana.shp")
bm <-  st_read("temp/BIEN2_ranges/2/Abutilon_hirtum.shp")
class(am)
am$geometry
plot(am[[1]][1])

j <-  st_intersects(am,bm)

aa=sf::as_Spatial(am)
bb=sf::as_Spatial(bm)
cc <- rbind(aa,bb)

data(africa)


bbb <- polys2comm(dat = cc, species = "species")
bbb

plot(bbb$poly_shp)

bbb$poly_shp@proj4string

plot(bbb$poly_shp)
gg <- st_as_sf(bbb$poly_shp)

ggplot()+
  geom_sf(data=wrld)+
  #geom_sf(data=bm,fill="red")+
  #geom_sf(data=am,fill="blue")+
  geom_sf(data=sp_grid)
  

amer_ras <- raster()
# Establecer el "extent" del raster para quedarnos justamente en las coordenadas de América del Sur
extent(amer_ras) <- c(-110,-29,-56,14)
res(amer_ras) <- 1
plot(amer_ras)
am_raster <- terra::rasterize(am, amer_ras)

plot(am_raster)
plot(wrld,add=T)


pa_matrix <- rasterToPoints(am_raster)
raster2comm(am_raster)

ggplot()+
  #geom_sf(data=wrld)+
  geom_sf(data=sp[[2]],fill="blue")



cbind(as.data.frame(pa_matrix),am$species)


#### 

cc <- PAM$comm_dat
object.size(cc)
s <- as.data.frame(rowSums(cc))
s <- cbind(s, rownames(cc))
colnames(s) <-  c("s","cell")


object.size(as.matrix(cc))

gg <- left_join(sp_grid,s,by="cell")
head(gg)

aa= st_as_sf(gg,coords = c(2,3))
aa


ggplot(data = aa)+
  geom_sf(aes(color=s),shape=15,size=0.8)+
  scale_color_viridis_c()


############################################################

aa=list()

ptm <- proc.time()
for (i in 1:50){
  temp <-  st_read(paste0("temp/BIEN2_ranges/1/",file.list[i]))
  aa[[i]] <-  sf::as_Spatial(temp)
}
proc.time() - ptm

aa <- aa[501:503]

ee <- do.call(rbind, aa)
ee <- cbind(ee, 1)
names(ee) <-  c("binomial","gid","presence")

ee@data <-  ee@data[,-2]

ee <- cbind(ee, 1,1)
names(ee) <-  c("binomial","presence","origin","seasonal")

plot(ee)

ff <- lets.presab(ee, xmn = -93.56732, xmx = -64.77008,
            ymn = -17.8574, ymx = 17.57342 ,cover=0,resol = 1)
ff2 <- lets.presab(ee, xmn = -106.5944, xmx = -29.33156,
                  ymn = -35.70736, ymx = 23.66393 ,cover=0,resol = 1)

ff$Presence_and_Absence_Matrix
str(ff$Presence_and_Absence_Matrix)
str(ff2$Presence_and_Absence_Matrix)

ff3 <- merge(ff$Presence_and_Absence_Matrix, 
      ff2$Presence_and_Absence_Matrix, by = c("Longitude(x)", "Latitude(y)"),
      all=T)
str(ff3)

ff$Richness_Raster

PAM <-  phyloregion::polys2comm(dat=ee, species="binomial",res=1)
PAM2 <-  phyloregion::polys2comm(dat=ee, species="binomial",res=1)

str(PAM$comm_dat)
str(PAM2$comm_dat)


########################################

j

raster_layer <- raster(extent(shplist[[78]]), res = 2)
bb <- rasterize(shplist[[78]], raster_layer)
plot(bb)


tochecklist <- c()
for(z in 1:length(df)){
  c <-df[z,]
  
  count <- 0
  
  for (i in 1:length(c@polygons[[1]]@Polygons)) {
    area <- c@polygons[[1]]@Polygons[[i]]@area
    if (area <  0.0082) {
      count <- count + 1
    }
    if (area <  0.0082) {
      tochecklist[z] <- df[z,]$binomial
    }
    
    print(count)
  }
}


z=78

a=lets.presab(df[1:78,], xmn = -180, xmx = 180, ymn = -90, ymx = 90,
              cover=0,resol = 1)
a$Presence_and_Absence_Matrix

plot(df[78,])
plot(df[1,],add=T)

a <- phyloregion::polys2comm(dat=df, species="species",res=1)

sp_grid <-  as.data.frame(coordinates(a$poly_shp))
sp_grid <-  cbind(rownames(sp_grid), sp_grid)
colnames(sp_grid) <- c("cell","X","Y")
sp_grid$cell  <- paste0("v",sp_grid$cell)

plot(a$poly_shp)






plot(df[98,])
extent(df[98,])
df[98,]@bbox[1:4] = round(df[98,]@bbox,2)

c <- df[98,]

c@bbox[1:4] <- round(c@bbox)


lets.presab(c,xmn = -180, xmx = 180, ymn = -90, ymx = 90,
            cover=0,resol = 1,count=T)
plot(c)
str(c)

amer_ras <- raster()

extent(amer_ras) <- c(-180,180,-90,90)
res(amer_ras) <- 1

am_raster <- terra::rasterize(c, amer_ras)
plot(am_raster)
terra::extract(am_raster,c)


### work with the PAMs
am <-  read_csv("temp/BIEN2_ranges/1/PAM1.csv")
bm <-  read_csv("temp/BIEN2_ranges/7/PAM6.csv")

bm[,"Acalypha_stachyura"]

colnames(am)[1:2] =c("X","Y")
bm = bm[,-1]

wrld
cc = full_join(am,bm)
cc <- cc %>% replace(is.na(.), 0)
ccc  <- as.data.frame(rowSums(cc[,-c(1:2)]))
colnames(ccc) <- "S"

ccc = cbind(cc[,1:2],ccc)

ccc=st_as_sf(ccc,coords = c(1,2))
st_crs(ccc)= st_crs(wrld)

ggplot()+
  geom_sf(data=wrld)+
  geom_sf(data=ccc,aes(color=S),shape=15,size=1)

vv=df[c(1,7),]

dd <-  phyloregion::polys2comm(dat=vv, species="species",res=1)



fdir <- system.file("NGAplants", package="phyloregion")
files <- file.path(fdir, dir(fdir))
ras <- raster2comm(files) # Note, this function generates
# a list of two objects
head(ras[[1]])



s <- readRDS(system.file("nigeria/nigeria.rds", package="phyloregion"))
sp <- random_species(100, species=5, shp=s)
pol <- polys2comm(dat = sp, species = "species")
head(pol[[1]])


s <- readRDS(system.file("nigeria/nigeria.rds", package = "phyloregion"))

set.seed(1)
m <- data.frame(sp::spsample(s, 10000, type = "nonaligned"))
names(m) <- c("lon", "lat")
species <- paste0("sp", sample(1:1000))
m$taxon <- sample(species, size = nrow(m), replace = TRUE)

pt <- points2comm(dat = m, mask = s, res = 0.5, lon = "lon", lat = "lat",
                  species = "taxon") # Note, this generates a list of two objects
head(pt[[1]])

makePAM(path="C:/Users/gabri/Dropbox/papers/projetos/Moulatlet - BIEN/temp/BIEN2_ranges/1",
        count=T)

am <-  st_read(paste0("temp/BIEN2_ranges/7/",ff[i]))
raster_grid <- raster(extent(am), res = c(0.1, 0.1))
grid_points <- as(raster_grid, "SpatialPoints")
presence_absence <- over(grid_points, am)

presence_absence_matrix <- as.data.frame(matrix(0, nrow = nrow(grid_points), ncol = nrow(am)))
colnames(presence_absence_matrix) <- spatial_data$your_attribute_column_name

for (i in 1:nrow(presence_absence)) {
  if (!is.na(presence_absence$your_attribute_column_name[i])) {
    presence_absence_matrix[i, presence_absence$your_attribute_column_name[i]] <- 1
  }
}
# Load the required libraries
library(terra)
library(sp)
library(raster)

# Load your SpatialPolygonsDataFrame (replace 'your_shapefile.shp' with your file)
spatial_data <- temp
plot(temp)

# Convert the spatial data to a terra object
terra_data <- vect(spatial_data)
plot(terra_data)

# Create a raster template with the same extent and resolution as the spatial data
template_raster <- rast(resolution = c(1, 1))


# Extract presence-absence information by rasterizing the spatial data
presence_absence <- rasterize(terra_data, template_raster, field = 1)
plot(presence_absence)
presence_absence[which(is.nan(terra::values(presence_absence)))] =0

aa = as.polygons(presence_absence)

one <- (presence_absence==1)
zero <- (presence_absence==0)
zero=as.polygons(zero)
plot(one)
one=as.polygons(one)
plot(zero)
str(presence_absence)

onepa =terra::extract(presence_absence,aa,xy=T)
zeropa =terra::extract(presence_absence,zero,xy=T)

vv =rbind(onepa,zeropa)
dim(onepa)
head(onepa)
plot(aa)

cc <- st_as_sf(onepa,coords = c(3,4))
vect(cc)
plot(cc)

wrld <- ne_countries(scale = 'small', returnclass = "sf")
st_crs(cc)= st_crs(wrld)

ggplot()+
  geom_sf(data=wrld)

save.image("workspace_27_09_23.RData")
