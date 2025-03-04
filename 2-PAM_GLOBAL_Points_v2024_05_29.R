# Load the required libraries
library(terra)         # CRAN v1.7-29
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0


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


# Load the object to be used as reference for the extraction

w <- vect("wrld_moll.shp")

# get the file names -----------------------------------------------------

file.list <- ptifs$Name

#split the data in several datasets

groups <- split(file.list, rep(1:100, each = length(file.list) /100, 
                                    length.out = length(file.list)))
namesg <- split(ptifs_sp, rep(1:100, each = length(file.list) /100, 
                               length.out = length(file.list)))


# #Loop1:  over each set of species ---------------------------------------

for (j in 1:length(groups)){ 
#for (j in 1:2){ #LOOP1
print(j)
  
glist <- groups[[j]] 
namel <- namesg[[j]] 
  

# Loop2: over each raster extract the values ----------

#create the result object
extList <- list()

for (i in 1:length(glist)){ #LOOP 2
#for (i in 1:100){ #LOOP 2

ras <- rast(unzip("temp/points_v1.zip",files = glist[i] ))
print(i)

#Extract presence-absence information by rasterizing the spatial data ---
    
names(ras)  <- namel[i]
values <- terra::extract(ras,w,touch=T,xy=T)
values <- values |> drop_na() |> pivot_longer(cols = namel[i])

#save the output
extList[[i]] <- values


} # end of loop 2
lf <- extList |> map(as_tibble) %>%
  reduce(bind_rows)

saveRDS(lf, paste0("products/PAMs_global_moll/","PAM_POINT","_",j,".rds"),compress = T)


# make the raster  
lfw <- lf |> 
  group_by(x,y) |> 
  summarise(S = sum(value),.groups = "keep")

  PAMrast <- as_spatraster(lfw, xycols=1:2)
  writeRaster(PAMrast , paste0("products/PAMrast_global_moll/","PAMraster_POINT","_",j,".tif"), overwrite=TRUE)

} #LOOP 1

