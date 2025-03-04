
############################################################
#                                                          #
#   Script to generate the simulated climatic layers       #
#                                                          #
############################################################

# script donwnloaded from Moore, T. E., Bagchi, R., Aiello‐Lammens, M. E., & Schlichting, C. D. (2018). 
# Spatial autocorrelation inflates niche breadth–range size relationships. 
# Global Ecology and Biogeography, 27(12), 1426-1436.
# https://onlinelibrary.wiley.com/doi/10.1111/geb.12818

# This script is not mine. It was made available with the above publication


# Load libraries ----------------------------------------------------------

library(geodata)
library(rnaturalearth)
library(geostatsp)
library(terra)
library(sf)
library(shapefiles)
library(INLA)
library(tidyverse)
library(usdm)
library(raster)
library(reshape2)
library(patchwork)
library(tidyterra)


# Import climate data -----------------------------------------------------

# importar el raster. Puede ser un raster con varios layers, uno
# para cada variable de interés. Yo hice con dos dimensiones solamentes  - PC1 y PC2 


mat.a <- rast("YOURRASTER.tif") 
res(mat.a)

if(F){
# check the spatial structure in the observed layer
v <- Variogram(mat.a)
plot(v)
}


# create matrix from the aggregated rasters
dat <- list()
for (i in 1:nlyr(mat.a)){
  temp <- terra::as.data.frame(mat.a[[i]],xy=T)
  colnames(temp) <-  c("x", "y", "z")
  dat[[i]] <- as.data.frame(temp)[with(as.data.frame(temp), order(y, decreasing = TRUE)),]
}


# Loop start --------------------------------------------------------------

for (j in 1:length(dat)){
  
# get the rows and columns into the correct order
  

xy_grid  <- reshape2::dcast(data =  dat[[j]], formula = y ~ x, value.var = "z") 
xy_grid  <- xy_grid[with(xy_grid, order(y, decreasing = TRUE)),]
y_vals   <- xy_grid[,1]
x_vals   <- as.numeric(colnames(xy_grid)[2:ncol(xy_grid)])
xy_grid  <- as.matrix(xy_grid)
xy_grid2 <- xy_grid[, -1]

mat.grid <- list(y_vals = y_vals, x_vals = x_vals, xy_grid2 = xy_grid2)

# Run the INLA model ------------------------------------------------------
# variables needed for the f function and the model

nrow = nrow(mat.grid$xy_grid2) 
ncol = ncol(mat.grid$xy_grid2) 
n = nrow*ncol
s.noise = 1
y = inla.matrix2vector(mat.grid$xy_grid2)
node = 1:n

formula = y ~ 1 + f(node, model="matern2d", nu=1, nrow=nrow, ncol=ncol,
                   hyper = list(range = list(param =c(1, 1),
                                             prior = "normal", param = c(0, 3)),
                                prec = list(param=c(1, 1)))) 


# If you need to check the initial values for optimization
inla.set.control.inla.default()

# get data into format for the model
mat.data = data.frame(y=y, node=node)


## fit the model
# NOTE: this may take some time to run if layer is large or at high resolution!

result <- inla(formula, family="gaussian", data=mat.data, verbose=TRUE,
            control.predictor = list(compute = TRUE),
            control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                            fixed = FALSE))),
            control.compute=list(cpo=FALSE), 
            control.inla = list(h = 1e-5, tolerance=1e-5),  
            keep=TRUE)

# Havard Rue suggested checking the mode.status to confirm that the model has been parameterized 
# appropriately. Ideally the mode.status should be 0
result$mode$mode.status 

# if this does not return 0, you can rerun the model using the previously estimated parameters as the starting points

if(F){
result = inla.rerun(result)
# and check the mode.status again
result$mode$mode.status
}

# store the parameter estimates to be used for simulation later

summ <- summary(result)

## plot the posterior mean for `predictor' and compare with the truth

# observed
#INLA:::inla.display.matrix(mat.grid$xy_grid2) 
#INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol))

# get predicted values into matrix
pred_mat <- inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)
colnames(pred_mat) <- mat.grid$x_vals
rownames(pred_mat) <- mat.grid$y_vals

# convert into a format that can be plotted in ggplot
pred_df <- as.data.frame(pred_mat)
colnames(pred_df)<- paste("V", colnames(pred_df), sep="_")

pred_df <- pred_df %>%
  rownames_to_column("y") %>%
  gather(x,z,starts_with("V"))

list_x <- unlist(strsplit(pred_df$x, split='_', fixed=TRUE))
pred_df$x <- as.numeric(list_x[c(seq(from=2, to=length(list_x), by = 2))])
pred_df$y <- as.numeric(pred_df$y)

#plot it
ggplot()+
  geom_raster(data = pred_df, aes(x = x, y=y, fill = z))

# get the observed data into the same set up 
obs_mat <- mat.grid$xy_grid2
colnames(obs_mat) <- mat.grid$x_vals
rownames(obs_mat) <- mat.grid$y_vals

obs_df <- as.data.frame(obs_mat)
colnames(obs_df)<- paste("V", colnames(obs_df), sep="_")

obs_df <- obs_df %>%
  rownames_to_column("y") %>%
  gather(x,z,starts_with("V"))

list_x_obs <- unlist(strsplit(obs_df$x, split='_', fixed=TRUE))
obs_df$x <- as.numeric(list_x_obs[c(seq(from=2, to=length(list_x_obs), by = 2))])
obs_df$y <- as.numeric(obs_df$y)


# plot it in ggplot
# observed
ggplot()+
  geom_raster(data = obs_df, aes(x = x, y=y, fill = z)) 

# remove the values that are NA in the observed dataset
pred_df2 <- pred_df[!is.na(obs_df[, 3]),]

# modeled
ggplot()+
  geom_raster(data = pred_df2, aes(x = x, y=y, fill = z)) 


# now plot the relationship between observed and predicted values
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
abline(b=1, a =0)
#summary(lm(inla.matrix2vector(pred_mat)~inla.matrix2vector(obs_mat)) )

if(F){
#-------------------------------------------------------------------------------------------------#
# this is to rerun the models with a different parameter value, i.e., setting
# nu = 2 and compare to 'result'

formula2= y ~ 1 + f(node, model="matern2d", nu=2, nrow=nrow, ncol=ncol,
                    hyper = list(range = list(param =c(1, 1),
                                              prior = "normal", param = c(0, 3)),
                                 prec = list(param=c(1, 1)))) 

result2=inla(formula2, family="gaussian", data=mat.data, verbose=TRUE,
             control.predictor = list(compute = TRUE),
             control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                             fixed = FALSE))),
             control.compute=list(cpo=FALSE), 
             control.inla = list(h = 1e-5, tolerance=1e-5),  
             keep=TRUE)

result2 = inla.rerun(result2)

summ2 <- summary(result2)
summ2

pred_mat2 <- inla.vector2matrix(result2$summary.linear.predictor$mean,nrow,ncol)
colnames(pred_mat2) <- mat.grid$x_vals
rownames(pred_mat2) <- mat.grid$y_vals


# now plot the relationship between the two models with different nu values
par(mfrow=c(1, 1))
#dev.new()
#dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(pred_mat2)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
}


# Simulating the climate data based on the INLA model ---------------------

#make a distance matrix
dmat <- as.matrix(dist(dat[[j]]))

# get parameter estimates from the INLA model
summ
# use these to set values below
beta0 <- -0.2157
sigma2e <- 1 / summ$hyperpar$mean[1]     # precision for gaussian observations
sigma2x <- 1 / summ$hyperpar$mean[2]     # precision for node
kappa <- sqrt(8) / summ$hyperpar$mean[3] # range for the node

# use value of nu specified for model
nu <- 1

mcor <- as.matrix(2 ^ (1 - nu) * (kappa * dmat) ^ nu *  besselK(dmat * kappa, nu) / gamma(nu))
diag(mcor) <- 1
mcov <- sigma2e * diag(ncol(dmat)) + sigma2x * mcor

# convert varcovar matrix to sd matrix
# NOTE: computationally intensive!!

L <- chol(mcor) # convert var covar to sd matrix

# sort the values of the real gradient in ascending order...
tdat <- data.frame(dat[[j]])
vals = as.numeric(sort(tdat$z))

#set up empty matrix and fill with predicted values
preds <- matrix(nrow = nrow(tdat), ncol = 50)
dim(preds)
for (i in 1:50) {
  preds[, i] <-
    beta0 + drop(rnorm(ncol(dmat)) %*% L) # simulate data
  print(i) # can suppress this
  
  # replace values in simulated layers with real ones in rank order
  # this is so that all layers have the same values in them 
  # (see Chapman 2010, citation in the paper)
  
  preds[, i] = vals[rank(preds[, i], ties.method = "random")]
}

preds <- as.data.frame(preds)
colnames(preds) <- paste("y", 1:50, sep = "")

# add coordinates
preds$x <- tdat$x
preds$y <- tdat$y

# check histograms to make sure the data match
hist(preds$y2)
hist(dat[[1]]$z, add = TRUE, col = "red")

if(F){
# example of observed vs expected
g1 <- ggplot() +
  geom_spatraster(data = mat.a , aes(fill= wc2.1_2.5m_bio_1)) +
  ggtitle("Observed") 
g1

# p
g2 <- ggplot(data = preds, aes(x = x, y = y, fill = y20)) +
  geom_raster() +
  ggtitle("Predicted")
g2


tt <- preds |> select(y25,y1000,x,y)
tt <- tt |> relocate(y25, .after=y) |> relocate(y1000,.after=y)
tttt <- rast(tt, type="xyz",crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

g2 <- ggplot( ) +
  geom_spatraster(data = tttt,aes(fill = y25)) +
  ggtitle("Predicted") +
  scale_fill_gradientn(colors = terrain.colors(10))
g2

g1+g2

max(dat[, 2])-min(dat[, 2])
}

# you can now write the scaled predicted layers to CSV
write.csv(preds, paste0("products/INLA/","predicted_inla_scaled_Americas_CV",j,".csv"))
}
save.image("INLAMERICA_CV.Rdata")


##%######################################################%##
#                                                          #
####                 Create tif layers                  ####
#                                                          #
##%######################################################%##
library(terra)
library(tidyverse)

files <- dir_ls(here("products","INLA", "predicted 50 km"))
pred1am <- read_csv(files[[1]])
pred2am <- read_csv(files[[2]])
pred1as <- read_csv(files[[3]])
pred2as <- read_csv(files[[4]])

m <- rast("Clim_scaled_world_moll_PCA.tif")


if(F){
# unscale values of MAP before saving and extracting for species
  
  preds.s <- preds
for (i in 1:1000){
  
  preds.s[, i] <- preds.s[, i]*sd(values(mat.a), na.rm = TRUE) + mean(values(mat.a), na.rm = TRUE)
}


# pred

# you can now write the unscaled predicted layers to CSV
write.csv(preds.s, "products/INLA/predicted_inla_unscaled.csv")

head(preds[, 1:10])

# check 
# you can use this method to make the predicted maps
pred.1 <- mat.a
values(pred.1)[!is.na(values(pred.1))] <- preds[, 100]
plot(pred.1) # first predicted layer
# compare
plot(mat.a) # observed layer
# we'll formalize this below
}

if(F){
# create raster layers from simulated data

sims <- list()
sims_grid <- list()
sims_r <-list()
i=1

for (i in 1:10){
  
  sims[[i]] <- as.data.frame(pred1[, c("x", "y", paste("y", i, sep = ""))])
  sims[[i]] <- sims[[i]][
    with(sims[[i]], order(y, decreasing = TRUE)),
  ]
  sims_grid[[i]] <-
    dcast( data =  sims[[i]], 
           formula = y ~ x, value.var = paste("y", i, sep = "")) 
  
  sims_grid[[i]] <- sims_grid[[i]][
    with(sims_grid[[i]], order(y, decreasing = TRUE)),
  ]
  
  sims_grid[[i]] <- as.matrix(sims_grid[[i]])
  sims_grid[[i]] <- sims_grid[[i]][, -1]
  sims_r[[i]] <- raster::raster(sims_grid[[i]], template = mat.a)
} 
}  





#### not in use ####

#sims_r now holds raster layers of each predicted env layer (should be 1000)

# example
plot(sims_r[[3]])

#check that details are the same (mins and maxes, etc)
mat.a
sims_r[[1]]



#--------------------------------------------------------------------------------------------------#

# generate completely randomized climate layers 

#--------------------------------------------------------------------------------------------------#
# create completely randomized layers of climate vars and save them as a csv. 

# 6858 cells in temp layer
# 6936 cells in rainfall layer

n.cell <- length(na.omit(values(mat.a)))


# first randomize data in 1000 maps
mat.a.r <- mat.a
random_maps <- list()

for (i in 1:10){
  
  random_maps[[i]] <- mat.a.r
  values(random_maps[[i]])[!is.na(values(random_maps[[i]]))] <- 
    values(random_maps[[i]])[!is.na(values(random_maps[[i]]))][sample(n.cell, n.cell)] # change numbers!!
}

# plot example
par(mfrow = c(2, 2))
plot(random_maps[[1]], main = "random")
plot(mat.a, main = "observed")
hist(random_maps[[1]], main = "random") # histograms should match
hist(mat.a, main = "observed") # histograms should match

# extract data for all pel species
pel.preds.rand <- as.data.frame(matrix(nrow = length(pel.sub$TAXNAME), ncol = 1000))
colnames(pel.preds.rand) <- paste("rep", 1:1000, sep = '.')

for (i in 1:1000) {
  
  pel.preds.rand[, i] <- raster::extract(random_maps[[i]], pel.sub[, c("LONGITUDE", "LATITUDE")])
  
  print(i)
}

pel.preds.rand$TAXNAME <- pel.sub$TAXNAME
write.csv(pel.preds.rand, "results/pel_pred_rand_scaled.csv")
# this will be used in observed vs simulated NB R code


