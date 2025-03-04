############################################################
#                                                          #
#  Script to generate 2d niches using observed  layers     #
#                                                          #
############################################################


# Load libraries ----------------------------------------------------------
library(here)
library(fs)
library(tidyverse)
library(hypervolume)
library(terra)
library(tidyterra)
library(future)
library(doFuture)
library(tictoc)
library(geometry)



###########################################################
#                                                         #
#            Start loop     RB                            #
#                                                         #
###########################################################

#load climatic data
rm <- rast("/home/u8/gmoulatlet/hypervolume/simulated/Clim_50km_scaled_moll.tif")
rmlist <- dir_ls("/home/u8/gmoulatlet/hypervolume/simulated/data",regexp = ".tif")

r <- list()
for (t in 1:length(rmlist)){
  v <- rast(rmlist[t])
  names(v) <- c("PC1","PC2")
  r[t] <- v
}
rm1 <- sprc(r)  
rm2 <- sprc(rm,rm1)

# Calculate niche breadth for each set of simulated layers -----------------
#load data
PAMs <- list.files("/home/u8/gmoulatlet/hypervolume/simulated/PAMS", "rds$")
PAMs_po <- PAMs[grep('RB',PAMs)]

tic()
# star the loop parallelization -------------------------------------------
for(j in 1:length(PAMs_po)) { #LOOP 1
#  for(j in 1:2) { #LOOP 1
  print(paste0("new PAM RB -",j))
  
  #select the PAM and extract the clim values
  PAM_plants <- as.data.frame(read_rds(paste0('/home/u8/gmoulatlet/hypervolume/simulated/PAMS/', PAMs_po[j])))
  PAM_plants <- PAM_plants |> rename(rangeID=ID) |> mutate(ID=seq(1:nrow(PAM_plants)))
  
  #get unique species
  sp <- PAM_plants |> dplyr::select(name) |> distinct()
  sp <- sp$name
  
  #prepare the output
  df <- data.frame(matrix(data=NA,nrow=300,ncol=52))
  colnames(df)[1:2] <- c('species','convex')
  for(u in 3:ncol(df)){
    colnames(df)[u] <- paste0("convexSimu",u-2)
  }
  
  for(k in 1:length(rm1)){
  #for(k in 1:2){
  clim_i <- terra::extract(rm2[k], PAM_plants[,c('x','y')], method = 'simple', na.rm = T,
                           ID=T,touch=T) 
  PAM_plants2 <- PAM_plants |> left_join(clim_i,by=c('ID'))
  
  sp <- PAM_plants2 |> dplyr::select(name) |> distinct()
  sp <- sp$name
  
 
  # loop over species -------------------------------------------------------
  for(i in 1:length(sp)) { #LOOP 1
  #for(i in 1:6) { #LOOP 1
    print(i)
    
    #prepare the species
    PAM_plants2 |> dplyr::filter(name==sp[i]) |> 
      dplyr::select(ID,x,y,name,PC1,PC2)-> temp;dim(temp)
    
    temp <- temp |> dplyr::select(name, PC1,PC2) |> distinct() ;dim(temp)
    
    df[i, 1] <- sp[i]
    
    #remove if the species has less than 20 occurrences - WONT USE THIS WITH THE CONVEX HULL, SO I AM ADDING <0
    if (length(temp[, 1]) < 3) { # loop removal species with low occurrences loop 3
      df[i, 1+k]=0
    } else{
      
      
      # extract from original data ----------------------------------------------
      clim_i0 <- temp |> dplyr::select(PC1,PC2)
      
      #delete NAs, if there are any
      #from orginal data
      clim_ina0 <- which(is.na(rowSums(clim_i0)))
      
      if (length(clim_ina0) > 0) {
        clim_i02 <- clim_i0[-clim_ina0, ]
      } else {
        clim_i02 = clim_i0
      }
      
      clim_i02  <- clim_i02 |> distinct()
      
      #make the climatic extracted values into a matrix
      m0 <- as.matrix(clim_i02)
      
      # calculate the hypervolme ------------------------------------------------
      #for the original data
      if (length(m0[, 1]) < 3) {
        df[i, 1+k] =0
      } else{
        
        gv0 <-  convhulln(m0, options="FA")
        df[i, 1+k] <- gv0$vol
      }
        
        #gv0 <- hypervolume_gaussian(m0, name = sp[i], verbose = FALSE)
        # df[i, 3]  <- get_volume(gv0)[1]
        
    }  
    } #close loop 3 -SPECIES
  } # close loop 2 - RASTERS
  #toc()
  write_rds(df, paste0('/home/u8/gmoulatlet/hypervolume/simulated/products/cvSimuRB_',j,'.rds'))  
    
} # close loop 1 - PAMS
toc()

############################################################
#                                                          #
#            Start loop     POINTS                         #
#                                                          #
############################################################


# Calculate hypervolumes for each set of simulated layers -----------------
#load data
PAMs <- list.files("/home/u8/gmoulatlet/hypervolume/simulated/PAMS", "rds$")
PAMs_po <- PAMs[grep("POINT",PAMs)]


tic()
# star the loop parallelization -------------------------------------------
for(j in 1:length(PAMs_po)) { #LOOP 1
  #  for(j in 1:2) { #LOOP 1
  print(paste0("new PAM POINTS-",j))
  
  #select the PAM and extract the clim values
  PAM_plants <- as.data.frame(read_rds(paste0('/home/u8/gmoulatlet/hypervolume/simulated/PAMS/', PAMs_po[j])))
  PAM_plants <- PAM_plants |> rename(rangeID=ID) |> mutate(ID=seq(1:nrow(PAM_plants)))
  
  #get unique species
  sp <- PAM_plants |> dplyr::select(name) |> distinct()
  sp <- sp$name
  
  #prepare the output
  df <- data.frame(matrix(data=NA,nrow=300,ncol=52))
  colnames(df)[1:2] <- c('species','convex')
  for(u in 3:ncol(df)){
    colnames(df)[u] <- paste0("convexSimu",u-2)
  }
  
  for(k in 1:length(rm1)){
    #for(k in 1:2){
    clim_i <- terra::extract(rm2[k], PAM_plants[,c('x','y')], method = 'simple', na.rm = T,
                             ID=T,touch=T) 
    PAM_plants2 <- PAM_plants |> left_join(clim_i,by=c('ID'))
    
    sp <- PAM_plants2 |> dplyr::select(name) |> distinct()
    sp <- sp$name
    
    
    # loop over species -------------------------------------------------------
    for(i in 1:length(sp)) { #LOOP 1
      #for(i in 1:6) { #LOOP 1
      print(i)
      
      #prepare the species
      PAM_plants2 |> dplyr::filter(name==sp[i]) |> 
        dplyr::select(ID,x,y,name,PC1,PC2)-> temp;dim(temp)
      
      temp <- temp |> dplyr::select(name, PC1,PC2) |> distinct() ;dim(temp)
      
      df[i, 1] <- sp[i]
      
      #remove if the species has less than 20 occurrences - WONT USE THIS WITH THE CONVEX HULL, SO I AM ADDING <0
      if (length(temp[, 1]) < 3) { # loop removal species with low occurrences loop 3
        df[i, 1+k]=0
      } else{
        
        
        # extract from original data ----------------------------------------------
        clim_i0 <- temp |> dplyr::select(PC1,PC2)
        
        #delete NAs, if there are any
        #from orginal data
        clim_ina0 <- which(is.na(rowSums(clim_i0)))
        
        if (length(clim_ina0) > 0) {
          clim_i02 <- clim_i0[-clim_ina0, ]
        } else {
          clim_i02 = clim_i0
        }
        
        clim_i02  <- clim_i02 |> distinct()
        
        #make the climatic extracted values into a matrix
        m0 <- as.matrix(clim_i02)
        
        # calculate the hypervolme ------------------------------------------------
        #for the original data
        if (length(m0[, 1]) < 3) {
          df[i, 1+k] =0
        } else{
          
          gv0 <-  convhulln(m0, options="FA")
          df[i, 1+k] <- gv0$vol
        }
        
        #gv0 <- hypervolume_gaussian(m0, name = sp[i], verbose = FALSE)
        # df[i, 3]  <- get_volume(gv0)[1]
        
      }  
    } #close loop 3 -SPECIES
  } # close loop 2 - RASTERS
  #toc()
  write_rds(df, paste0('/home/u8/gmoulatlet/hypervolume/simulated/products/cvSimuPoints_',j,'.rds'))  
  
} # close loop 1 - PAMS
toc()



############################################################
#                                                          #
#            Start loop     PPM                            #
#                                                          #
############################################################


# Calculate hypervolumes for each set of simulated layers -----------------
#load data
PAMs <- list.files("/home/u8/gmoulatlet/hypervolume/simulated/PAMS", "rds$")
PAMs_po <- PAMs[grep("PPM",PAMs)]


tic()
# star the loop parallelization -------------------------------------------
for(j in 1:length(PAMs_po)) { #LOOP 1
  #  for(j in 1:2) { #LOOP 1
  print(paste0("new PAM PPM -",j))
  
  #select the PAM and extract the clim values
  PAM_plants <- as.data.frame(read_rds(paste0('/home/u8/gmoulatlet/hypervolume/simulated/PAMS/', PAMs_po[j])))
  PAM_plants <- PAM_plants |> rename(rangeID=ID) |> mutate(ID=seq(1:nrow(PAM_plants)))
  
  #get unique species
  sp <- PAM_plants |> dplyr::select(name) |> distinct()
  sp <- sp$name
  
  #prepare the output
  df <- data.frame(matrix(data=NA,nrow=300,ncol=52))
  colnames(df)[1:2] <- c('species','convex')
  for(u in 3:ncol(df)){
    colnames(df)[u] <- paste0("convexSimu",u-2)
  }
  
  for(k in 1:length(rm1)){
    #for(k in 1:2){
    clim_i <- terra::extract(rm2[k], PAM_plants[,c('x','y')], method = 'simple', na.rm = T,
                             ID=T,touch=T) 
    PAM_plants2 <- PAM_plants |> left_join(clim_i,by=c('ID'))
    
    sp <- PAM_plants2 |> dplyr::select(name) |> distinct()
    sp <- sp$name
    
    
    # loop over species -------------------------------------------------------
    for(i in 1:length(sp)) { #LOOP 1
      #for(i in 1:6) { #LOOP 1
      print(i)
      
      #prepare the species
      PAM_plants2 |> dplyr::filter(name==sp[i]) |> 
        dplyr::select(ID,x,y,name,PC1,PC2)-> temp;dim(temp)
      
      temp <- temp |> dplyr::select(name, PC1,PC2) |> distinct() ;dim(temp)
      
      df[i, 1] <- sp[i]
      
      #remove if the species has less than 20 occurrences - WONT USE THIS WITH THE CONVEX HULL, SO I AM ADDING <0
      if (length(temp[, 1]) < 3) { # loop removal species with low occurrences loop 3
        df[i, 1+k]=0
      } else{
        
        
        # extract from original data ----------------------------------------------
        clim_i0 <- temp |> dplyr::select(PC1,PC2)
        
        #delete NAs, if there are any
        #from orginal data
        clim_ina0 <- which(is.na(rowSums(clim_i0)))
        
        if (length(clim_ina0) > 0) {
          clim_i02 <- clim_i0[-clim_ina0, ]
        } else {
          clim_i02 = clim_i0
        }
        
        clim_i02  <- clim_i02 |> distinct()
        
        #make the climatic extracted values into a matrix
        m0 <- as.matrix(clim_i02)
        
        # calculate the hypervolme ------------------------------------------------
        #for the original data
        if (length(m0[, 1]) < 3) {
          df[i, 1+k] =0
        } else{
          
          gv0 <-  convhulln(m0, options="FA")
          df[i, 1+k] <- gv0$vol
        }
        
        #gv0 <- hypervolume_gaussian(m0, name = sp[i], verbose = FALSE)
        # df[i, 3]  <- get_volume(gv0)[1]
        
      }  
    } #close loop 3 -SPECIES
  } # close loop 2 - RASTERS
  #toc()
  write_rds(df, paste0('/home/u8/gmoulatlet/hypervolume/simulated/products/cvSimuPPM_',j,'.rds'))  
  
} # close loop 1 - PAMS
toc()
