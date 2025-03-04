# Load packages -----------------------------------------------------------
library(tidyverse)
library(arrow)
library(tictoc)
library(here)
library(fs)
library(terra)

if(F){
  

# Load files --------------------------------------------------------------
pamlist <- list.files(path = paste0("products/PAMs_global_moll/"),pattern = c("*.rds"), full.names = TRUE)
length(pamlist)


# read and process points only --------------------------------------------
plist <- pamlist[grep("POINT",pamlist)]

#load metadata points
mdp <- read_rds(here("products",
                     "areas_midpoints_global_moll",
                     "area_midpoints_POINTS_regionsGLOBAL.rds"))
mdp <- mdp |>
  rename(name=species) |> 
  select(name,region,band)

data_list <- map(plist,read_rds)
data_list <- map_df(data_list,~ mutate(.x, type = "point"))
data_list <- dplyr::left_join(data_list,mdp,by="name")

write_parquet(data_list, "products/PAMs_global_moll_sum/POINTS.parquet")

# read and process ppm only -----------------------------------------------
plist <- pamlist[grep("PPM",pamlist)]

num_parts <- 50
groups <- cut(seq_along(plist), breaks = num_parts, labels = FALSE)
split_vector <- split(plist, groups)

#load metadata ppm
mdppm <- read_rds(here("products",
                     "areas_midpoints_global_moll",
                     "area_midpoints_PPM_regionsGLOBAL.rds"))
mdppm <- mdppm |>
  rename(name=species) |> 
  select(name,region,band)


plan(multisession)
foreach(j=41:50, .combine = c) %dofuture%{
  data_list <- purrr::map(split_vector[[j]],read_rds)
  data_list <- map_df(data_list,~ mutate(.x, type = "ppm"))
  data_list <- dplyr::left_join(data_list,mdppm,by="name")
  write_parquet(data_list, 
                paste0("products/PAMs_global_moll_sum/PPM",j,".parquet")
                )
  rm(data_list)
}


# read and process rb only -----------------------------------------------
#load metadata rb
alist <- list.files(path = here("products","areas_midpoints_global_moll"),pattern = c("*.rds"), full.names = TRUE)
alistrb <- alist[grep("RANGEBAGGING",alist)]
mdrb <- map_df(alistrb,read_rds)

plist <- pamlist[grep("RB",pamlist)]

num_parts <- 50
groups <- cut(seq_along(plist), breaks = num_parts, labels = FALSE)
split_vector <- split(plist, groups)

mdrb <- mdrb |>
  rename(name=species) |> 
  select(name,region,band)

plan(multisession)
foreach(j=1:50, .combine = c) %dofuture%{
  data_list <- purrr::map(split_vector[[j]],read_rds)
  data_list <- map_df(data_list,~ mutate(.x, type = "rb"))
  data_list <- dplyr::left_join(data_list,mdrb,by="name")
  write_parquet(data_list, 
                paste0("products/PAMs_global_moll_sum/RB",j,".parquet")
  )
  rm(data_list)
}

}
# Create the dataset ------------------------------------------------------
# Load taxonomic database -------------------------------------------------
taxo <- read_csv(here("temp","range_model_species_20230524_all.csv"))

taxo[which(taxo$higher_plant_group=="gymnosperms (non-conifer)"),"higher_plant_group"] <- 
  "gymnosperms"
taxo[which(taxo$higher_plant_group=="gymnosperms (conifers)"),"higher_plant_group"] <- 
  "gymnosperms"
taxo <- taxo |> select(species_nospace,higher_plant_group)
taxo <- taxo |> rename(name=species_nospace)

r <-
  arrow::open_dataset(
    sources = here("products",
                   "PAMs_global_moll_sum"),
    format = "parquet")


r2 <- r |> 
  dplyr::select(-value) |>  
  #dplyr::filter(type=="point") |> 
  left_join(taxo, by="name") |> 
  dplyr::group_by(type,band,region) |> 
  write_dataset(here("products",
                     "PAMs_global_moll_sum",
                     "ranges_ds"),format = "parquet")
  


# Load the recent created dataset ------------------------------------------------------------
ds <- open_dataset("products/PAMs_global_moll_sum/ranges_ds/", format = "parquet")



dim(ds)

dir_tree(here("products","PAMs_global_moll_sum","ranges_ds"))
total_size <- sum(dir_info(
  here("products","PAMs_global_moll_sum","ranges_ds"), 
  recurse = TRUE)$size, na.rm = TRUE); total_size

if(F){
#Associate elevations -------------------------------------------------------------------
  #this code has to run on HPC
  ele <- rast("Elevation_world_moll.tif")
  names(ele) <- "elevation"
  
  y <- 2
  i = y
  #load files 
  
  parquetlist <- list.files(path = paste0("products/PAMs_global_moll_sum/"),pattern = c("*.parquet"), full.names = TRUE)
  
  a <- read_parquet(parquetlist[i])
  a |> select(name) |> distinct() -> sps
  
  tic()
  resu <- data.frame(name=NA, min=NA, max=NA, mean=NA)
  for (j in 1:length(sps$name)){
    #for (j in 1:5){
    print(j)
    s <- sps$name[j]
    print(s)
    temp <-  a |> select(x,y,name) |> filter(name==s) |> mutate(IDele=seq_len(n()))
    ele_i <- terra::extract(ele, temp[,c("x","y")], method = "bilinear", na.rm = T,
                            ID=T)  |> rename(IDele=ID)
    t <- temp |> left_join(ele_i,by="IDele") |> 
      mutate(min=min(elevation),max=max(elevation),mean=mean(elevation)) |> 
      select(name,min,max,mean) |> 
      distinct()
    resu[j,] <- t
  }
}
############################################################################