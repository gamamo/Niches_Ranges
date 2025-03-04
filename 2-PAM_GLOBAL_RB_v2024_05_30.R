# Load the required libraries
library(terra)         # CRAN v1.7-29
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0

#process points #####
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

# Create a raster template with the same extent and resolution as the spatial data

w <- vect("wrld_moll.shp")

# get the file names -----------------------------------------------------

file.list <- rbtifs$Name
groups <- split(file.list, rep(1:500, each = length(file.list) /500, 
                               length.out = length(file.list)))
namesg <- split(rbtifs_sp, rep(1:500, each = length(file.list) /500, 
                              length.out = length(file.list)))

for (j in 473:length(groups)){ #LOOP1
print(paste0("new group"," = ", j))
  
glist <- groups[[j]] 
namel <- namesg[[j]] 
  
# Loop2: over each raster extract the values ----------

#create the result object
extList <- list()

for (i in 1:length(glist)){ #LOOP 2
 #for (i in 1:100){ #LOOP 2
    
ras <- rast(unzip("temp/rangebag_v1.zip",files = glist[i] ))
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

saveRDS(lf, paste0("products/PAMs_global_moll/","PAM_RB","_",j,".rds"),compress = T)

# make the raster  
lfw <- lf |> 
  group_by(x,y) |> 
  summarise(S = sum(value),.groups = "keep")

PAMrast <- as_spatraster(lfw, xycols=1:2)
writeRaster(PAMrast , paste0("products/PAMrast_global_moll/","PAMraster_RB","_",j,".tif"), overwrite=TRUE)
  
} #LOOP 1



