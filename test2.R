# Load the required libraries
library(terra)         # CRAN v1.7-29
library(raster)        # CRAN v3.6-20
library(sf)            # CRAN v1.0-13
library(tidyverse)     # CRAN v2.0.0
library(tidyterra)     # CRAN v0.4.0


# Create a raster template with the same extent and resolution as the spatial data
template_raster <- rast(resolution = c(1, 1),xmin=-180, xmax=180, 
                        ymin=-90, ymax=90)
coord <- crds(template_raster)

# get the file names -----------------------------------------------------

file.list <- list.files(path = paste0("temp/BIEN2_ranges/"),pattern = c("*.shp"),recursive = T)

# Group the elements into a list with groups with 989 elements each
grouped_list <- split(file.list, rep(1:ceiling(length(file.list)/2), 
                                     each = 989, length.out = length(file.list)))


# Loop1: over these groups -----------------------------------------------

for (j in 1:length(grouped_list)){
  
  g <- grouped_list[[j]]

  
# Create the PAM matrix --------------------------------------------------

community_df <- as.data.frame(matrix(NA, nrow=nrow(crds(template_raster))))
community_df[,1] <- crds(template_raster)[,1] # get the X coordinate
community_df[,2] <- crds(template_raster)[,2] # get the Y coordinate
colnames(community_df) <- c("X","Y")  
  
  
# Loop2: over each shapefile, rasterize, and extract the values ----------
  
for (i in 1:length(g)){

    shp<-  st_read(paste0("temp/BIEN2_ranges/",g[i]))
  
# Extract presence-absence information by rasterizing the spatial data ---
    
presence_absence <- terra::rasterize(shp, template_raster, field = 1,background=0,
                                     touches=T,cover=0)
names(presence_absence)  <- shp$species

values <- terra::extract(presence_absence,coord)

# Paste values in a community list ---------------------------------------

community_df[,(i+2)] <- cbind(values)
colnames(community_df)[(i+2)] <- names(values) 

} # end of loop 2


# Prepare the community table to export ----------------------------------

temp <- community_df[,-c(1:2)]

r_todelete <- which(rowSums(temp)==0)
PAM <- community_df[-r_todelete,]

saveRDS(PAM, paste0("products/PAM",j,".rds"))

# Make the richness table--------------------------------------------------

PAMrich <- as.data.frame(rowSums(community_df[,-c(1:2)]))
colnames(PAMrich ) <- "S"
PAMrich  = cbind(community_df[,c(1,2)],PAMrich )

# Make and export the raster ----------------------------------------------

PAMrast <- as_spatraster(PAMrich,xycols=1:2)
writeRaster(PAMrast , paste0("products/PAMrast",j,".tif"), overwrite=TRUE)


} # end of the loop 1





