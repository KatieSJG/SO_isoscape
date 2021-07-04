library(raster)
library(INLA)
library(maptools)
library(ggplot2)
library(cowplot)
library(rgeos)
library(SpatialEpi)
library(splancs) # for inout function

require(PBSmapping)
library(maps)
library(mapdata)
library(maptools)
library(sf)

polar = CRS("+init=epsg:3031")
longlat = CRS("+init=epsg:4326")

# 1. DATA INPUT ----
# load isotope data and get environmental data for the models
SO<-read.csv("POM data.csv",header = TRUE)
#Add in Season Columns
Jan<-subset(SO,SO$Month == "Jan")
Jan$Season<-1
Feb<-subset(SO,SO$Month == "Feb")
Feb$Season<-1
Mar<-subset(SO,SO$Month == "Mar")
Mar$Season<-2
Apr<-subset(SO,SO$Month == "Apr")
Apr$Season<-2
May<-subset(SO,SO$Month == "May")
May$Season<-3
Jun<-subset(SO,SO$Month == "Jun")
Jun$Season<-3
Jul<-subset(SO,SO$Month == "Jul")
Jul$Season<-3
Aug<-subset(SO,SO$Month == "Aug")
Aug$Season<-3
Sep<-subset(SO,SO$Month == "Sep")
Sep$Season<-3
Oct<-subset(SO,SO$Month == "Oct")
Oct$Season<-3
Nov<-subset(SO,SO$Month == "Nov")
Nov$Season<-4
Dec<-subset(SO,SO$Month == "Dec")
Dec$Season<-4

SO<-rbind(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)

#remove columns I dont want - to make a smaller data frame
SO$Source<-NULL
SO$mass..cm.<-NULL
SO$Length.cm.<-NULL
SO$C.N<-NULL
SO$Species<-NULL
SO$dN15<-NULL
#SO$dC13<-NULL
SO<-na.omit(SO)

range(SO$dC13)
SO<-subset(SO,SO$dC13 > -47)

# get points for sampling locations
x <- SO$Long
y <- SO$Lat
SO_xy <- cbind(x,y)
SO_pts <- SpatialPoints(SO_xy)


#Mean environ data ----

SST_ras<-raster("SST data.grd")
MLD_ras<-raster("MLD data.grd")
CHL_ras<-raster("CHL data.grd")
PPv_ras<-raster("PPV data.grd")
Dist_ras<-raster("Distance to coast data.grd")

#get rid of two PPV outliers 
values(PPv_ras)[values(PPv_ras) > 2000]=NA

#Worked out how to change resolution to 0.1 by 0.1 - do this in future:
#Make an empty raster so that resolution is 0.1 by 0.1 - use to resample env var rasters
blank <- raster(ncol=180, nrow=40, xmn=-180, xmx=180, ymn=-79.5, ymx=-39.5, crs = '+proj=longlat +ellps=sphere +a=1 +b=1')
# 
SST <- resample(SST_ras,blank)
CHL <- resample(CHL_ras,blank)
MLD <- resample(MLD_ras,blank)
PPv <- resample(PPv_ras,blank)
Dist <- resample(Dist_ras,blank)

# put the rasters in a stack and name them - makes extracting data easier
env_data <- stack(SST, CHL, MLD, PPv, Dist)
env_data <- scale(env_data) # standardising the values here so that they are standardised over the same mean and sd for fit and prediction
names(env_data) <- c("sst", "chl", "mld", "ppv", "dist")

# extract data for each sample location
env_data_pts <- extract(env_data, SO_pts)

# merge to the main data - this will just merge the standardised data onto the main dataframe
SO <- cbind(SO, env_data_pts)
SO <- na.omit(SO)

x <- SO$Long
y <- SO$Lat
SO_xy <- cbind(x,y)
SO_pts <- SpatialPoints(SO_xy)

# 2. CREATE THE MESH ----
# create the boundary for use in the mesh
#XY_km <- as.matrix(latlong2grid(cbind(SO$Long, SO$Lat)))

data("wrld_simpl")

land<- subset(wrld_simpl,NAME %in% c("Antarctica","Argentina","Australia","New Zealand"))
land <- gSimplify(land, tol=0.8)
land<-aggregate(land)
crs(land)<-"+init=epsg:4326"
land = spTransform(land, polar)
plot(land)

# 

#Try a buffer 

coords <- matrix(c(-1, -89,
                   -1, -90, 
                   1, -90,                          
                   1, -89,
                   -1, -89), 
                 ncol=2, byrow=TRUE)

P1 <- Polygon(coords)
Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")))
plot(Ps1)
crs(Ps1)<-"+init=epsg:4326"
ps1.polar = spTransform(Ps1, polar)
# 
ps1.polar2 <- gBuffer(ps1.polar, width = 6000000, byid=FALSE)
plot(ps1.polar2, col = "red")
plot(ps1.polar, col = 'black', add = TRUE)
plot(land, add = TRUE)

sea <- gDifference(ps1.polar2, land)
plot(sea)



# create the mesh for use in the spatial models

mesh_bound <- inla.mesh.2d(
  max.edge = (4*3.14)/180,
  boundary = as.inla.mesh.segment(sea),
  crs = inla.CRS("+proj=geocent +a=6370997 +b=6370997 +units=m +no_defs"))

crs(SO_pts)<-"+init=epsg:4326"
SO_pts2 <- sf::st_transform(sf::st_as_sf(SO_pts), "+proj=geocent +a=6370997 +b=6370997 +units=m +no_defs")
SO_pts3 <- as(SO_pts2, "Spatial")
plot(mesh_bound, asp = 1)
points(SO_pts3, pch=16,cex=0.5,col="red")
SO_pts4<-as.data.frame(SO_pts3)
SO_pts5<-SO_pts4[,1:2]


# 3. GENERATE PREDICTION DATA ----
# now we want to create the data to predict into - using the methods from https://www.stat.washington.edu/peter/591/INLA.html


# the aggregate line adjusts the resolution 
env_data_coarse <- aggregate(env_data, fact=1) 
plot(env_data_coarse, useRaster=T) # still reasonably fine
env_data_pred <- na.omit(as.data.frame(env_data_coarse, xy = TRUE))


range(SO$sst)
range(SO$chl)
range(SO$mld)
range(SO$ppv)
range(SO$dist)

pred_locs <- cbind(env_data_pred$x, env_data_pred$y)
pred_locs2 <- SpatialPoints(pred_locs)
crs(pred_locs2)<-"+init=epsg:4326"

pred_locs3 <- sf::st_transform(sf::st_as_sf(pred_locs2), "+proj=geocent +a=6370997 +b=6370997 +units=m +no_defs")
pred_locs4<- as(pred_locs3, "Spatial")
pred_locs5<-as.data.frame(pred_locs4)
pred_locs6<-as.matrix(pred_locs5[,1:2])

# now link the predicted locations to the mesh (and limit to the boundary we set in the mesh section)

#Don't need to limit to a boundary as the environmental data are already the correct size we want to project to

pred_grid <- inla.mesh.projector(mesh_bound, loc = pred_locs6)

# 
# 
# # limit the locations and environmental data to those points within the boundary (denoted by xy_in)
#pred_locs6 <- pred_locs6[xy_in,]
colnames(pred_locs6) <- c("x", "y")
#env_data_pred <- env_data_pred[xy_in,]

pred_pts <- SpatialPoints(pred_locs6)
plot(pred_pts) # check we are plotting into sensible places


# 4. ASSOCIATE OBSERVATION AND PREDICTION DATA LOCATIONS WITH MESH VERTICES (A MATRIX)
XY_km <- as.matrix(SO_pts5)


A_b<-inla.spde.make.A(mesh_bound, loc = XY_km)
spde_b<-inla.spde2.matern(mesh_bound, alpha = 2)
wb.index <- inla.spde.make.index('w', n.spde = spde_b$n.spde)

# we need to calculate the A for the predictions too

A_pred <- pred_grid$proj$A

# sample size
N_obs <- nrow(SO)
N_pred <- nrow(pred_locs6)
# Now we've got all the bits to go into the models, we can do the models:

# 4. MODELS ----

# 1. Create the stack

c_stack_fit <- inla.stack(
  tag = "Fit",
  data = list(y = SO$dC13),  
  A = list(1, 1, 1,1, 1,A_b),                  
  effects = list(   
    Intercept = rep(1, N_obs),
    X = SO[,c("sst", "mld","ppv", "dist")],
    Season = SO$Season,
    Year = SO$Year,
    StudyID = SO$StudyID,
    w = wb.index))

stack_pred <- inla.stack(data = list(y = NA), 
                         A = list(1, 1, A_pred), tag = "Predicted", 
                         effects = list(
                           Intercept = rep(1, N_pred), 
                           X = env_data_pred[,c("sst", "mld", "ppv", "dist")],
                           w = wb.index
                         ))

c_stack_all <- inla.stack(c_stack_fit, stack_pred)

# 3. Spatial model

#interactions: 

c_spatial_f <- y ~ -1 + Intercept + sst + mld + dist + ppv + sst:ppv + mld:ppv + dist:ppv +  f(Season, model = "iid") +   f(Year, model = "iid") + f(StudyID, model = "iid") + f(w, model = spde_b)

#no interactions: 

c_spatial_f <- y ~ -1 + Intercept + sst + mld + ppv + f(Season, model = "iid") +   f(Year, model = "iid") +  f(StudyID, model = "iid") + f(w, model = spde_b)

c_spatial_m <- inla(c_spatial_f, 
                    family = "gaussian", 
                    data = inla.stack.data(c_stack_fit),
                    control.compute = list(dic = TRUE), 
                    control.predictor = list(A = inla.stack.A(c_stack_fit)))


#3.a. model outputs
length(SO$dC13)

summary(c_spatial_m)
Fit1 <- c_spatial_m$summary.fitted.values$mean[1:3339]
E1   <- SO$dC13 - Fit1

plot(Fit1,SO$dC13)
cor.test(Fit1,SO$dC13)

#Homogeneity
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Fit1, y = E1)
abline(h = 0, v = 0)

#Normality
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(E1, breaks = 25)

#Independence due to model misfit
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = SO$Season, 
     y = E1)
abline(h = 0)

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = SO$Year, 
     y = E1)
abline(h = 0)

c_spatial_m$summary.random

c_spatial_m$summary.hyperpar

#Plot Random MArginals for C - Season
plot(c_spatial_m$marginals.random$Season$index.1, col = "grey", type = "l", lwd = 4, ylim=c(0,4), xlim=c(-3,3))
points(c_spatial_m$marginals.random$Season$index.2, col = "light blue", type = "l", lwd = 4)
points(c_spatial_m$marginals.random$Season$index.3, col = "goldenrod1", type = "l", lwd = 4)
points(c_spatial_m$marginals.random$Season$index.4, col = "darkseagreen2", type = "l", lwd = 4)

# 4. Predict
c_spatial_pred <- inla(c_spatial_f, 
                       family = "gaussian", 
                       data = inla.stack.data(c_stack_all),
                       control.compute = list(dic = TRUE), 
                       control.predictor = list(A = inla.stack.A(c_stack_all)),  
                       quantiles = NULL, 
                       control.results = list(return.marginals.random = FALSE, return.marginals.predictor = FALSE), 
                       control.inla = list(int.strategy = "eb"))

# 5. Extract the indices of the prediction nodes and extract mean and sd of the response 
id_pred <- inla.stack.index(c_stack_all, "Predicted")$data
c_mean <- c_spatial_pred$summary.fitted.values$mean[id_pred]
c_sd <- c_spatial_pred$summary.fitted.values$sd[id_pred]

#pred_locs2 <- pred_locs2[xy_in,]
pred_locs_final<-as.data.frame(pred_locs2)
colnames(pred_locs_final) <- c("x", "y")

c_res <- data.frame(pred_locs_final, c_mean, c_var = c_sd^2)

# 7. Plot the results - redo so plot as a raster 
c_mean_r<-subset(c_res, select = c("x", "y", "c_mean"))
coordinates(c_mean_r) <- ~x+y 
gridded(c_mean_r) <- TRUE
c_mean_raster<-raster(c_mean_r)


colfunc<-colorRampPalette(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695"))
colfunc(100)
plot(c_mean_raster,col=rev(colfunc(100)))
points(SO_xy,pch=16,cex=0.2,col="black")
map('worldHires',xlim=c(-180,180),ylim=c(-70,-40), col="white", fill=TRUE, add =TRUE)


c_var_r<-subset(c_res, select = c("x", "y", "c_var"))
coordinates(c_var_r) <- ~x+y 
gridded(c_var_r) <- TRUE
c_var_raster<-raster(c_var_r)

plot(c_var_raster,col=rev(colfunc(100)))
points(SO_xy,pch=16,cex=0.2,col="black")
map('worldHires',xlim=c(-13,9),ylim=c(48,62), col="white", fill=TRUE, add =TRUE)


#Save all rasters----
setwd(" folder ")
writeRaster(c_mean_raster, "c_mean_interactions",overwrite=TRUE)
writeRaster(c_var_raster, "c_var_interactions",overwrite=TRUE)

#Re run with no interaction model 
# Repeat with Nitrogen 


