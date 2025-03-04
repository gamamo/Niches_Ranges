# Load the required libraries
library(terra)         # CRAN v1.7-29
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0


#process points #####
ppm <- unzip("temp/ppm_v1/BinaryMaps.zip",list=T)

ppmt <- grep('\\.tif'  ,ppm$Name)
ppmtifs <- ppm[ppmt, ]
tp05   <- grep("/TP05__",ppmtifs$Name)
ppmtifs <- ppmtifs[tp05, ]
ppmocc <- grep("_occEco.tif",ppmtifs$Name)
ppmtifs <- ppmtifs[ppmocc, ]

# get only species names  
ppmtifs_sp <- gsub(".*Maps/","",ppmtifs$Name) 
ppmtifs_sp <- gsub("/TP05.*","",ppmtifs_sp)  
ppmtifs_sp <- data.frame(names=ppmtifs_sp)
ppmtifs_sp <- ppmtifs_sp$names


# Create a raster template with the same extent and resolution as the spatial data

w <- vect("wrld_moll.shp")

# get the file names -----------------------------------------------------

file.list <- ppmtifs$Name
groups <- split(file.list, rep(1:500, each = length(file.list) /500, 
                               length.out = length(file.list)))
namesg <- split(ppmtifs_sp, rep(1:500, each = length(file.list) /500, 
                               length.out = length(file.list)))

for (j in 1:length(groups)){ #LOOP1
  print(paste0("new group"," = ", j))
  
  glist <- groups[[j]] 
  namel <- namesg[[j]] 
  
# Loop2: over each raster extract the values ----------
  
#create the result object
extList <- list()
  
#for (i in 1:length(glist)){ #LOOP 2
  for (i in 1:50){ #LOOP 2

ras <- rast(unzip("temp/ppm_v1/BinaryMaps.zip",files = glist[i] ))
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
  
saveRDS(lf, paste0("products/PAMs_global_moll/","PAM_PPM","_",j,".rds"),compress = T)


# make the raster  
lfw <- lf |> 
  group_by(x,y) |> 
  summarise(S = sum(value),.groups = "keep")

PAMrast <- as_spatraster(lfw, xycols=1:2)
writeRaster(PAMrast , paste0("products/PAMrast_global_moll/","PAMraster_PPM","_",j,".tif"), overwrite=TRUE)
  
} #LOOP 1


