####	####	####	####
### CODE for "Smaller climatic niche shifts in invasive compared to non-invasive alien species"
## Olivia Bates, S?bastien Ollier, Cleo Bertelsmeier
# 06/2020
####	####	####	####


#### Packages required
library(ade4)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(dismo)
library(ecospat)
library(geosphere)
library(countrycode)
library(hypervolume)
library(outliers)
library(class)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggpubr)
library(grDevices)
library(extrafont)

#### Load data
rm(list = ls())
load("data_invasiveants.RData")

### Antmaps data
	# antmaps_bentities : names of the 546 bentities
	# antmaps_distri : information on the occurences of the 82 ant species
	# antmaps_iso3 : ISO3 code of the 546 bentities
	# antmaps_occ : occurences of the 82 species
	# antmaps_spol : spatial features of 428 bentities (no island)
	# antmaps_spol_i : spatial features of 107 bentities (islands)

### Bertelsmeier at al. (2017) data  
	# groups241 : dispersion groups of the 241 exotic species (1 : local ; 2 : regional ; 3 : trancontinental ; 4 : global)
	# species241 : names of the 241 exotic species
	# xy_map : rao spatial diversity and richness of the 13104 ant species
	# map : world map of the 244 countries

### Paper data
	# species82 : names of the 82 non-native species
	# groups82 : dispersion groups of the 82 non-native species
	# status82 : status (invasion or alien) for the 82 species
	# niche_results : niche shift statistics estimated in this paper (
	#--> used to built Table S1
	# sample_size : number of occurences for the 82 species
	#--> used to built Table S1

### Add climatic data
#--> 
path <- "C:/Users/Ollier/Desktop/CleoGroup/Data/Geomatic/Climatic/worldclim"
r <- raster::getData("worldclim", var = "bio", res = 2.5, path = path)
r <- r[[c(1, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)]]


#### Add local functions 
### For equivalency test (adapted from ecospat)
"overlap.eq.gen" <- function(repi, z1, z2) {
  sim.o.D <- NaN
  while (is.nan(sim.o.D)) {
    
    # overlap on one axis
    ##z1$sp is the environmental values for occurences of the species 
    occ.pool <- c(z1$sp, z2$sp)  # pool of random occurrences
    occ.pool <- occ.pool[sample(length(occ.pool))]
    row.names(occ.pool)<-c() # remove rownames 
    rand.row <- base::sample(1:length(occ.pool), length(z1$sp)) # random reallocation of occurrences to dataset
    sp1.sim <- occ.pool[rand.row]
    sp2.sim <- occ.pool[-rand.row]
    
    z1.sim <- ecospat.grid.clim.dyn(z1$glob, z1$glob1, data.frame(sp1.sim), R = length(z1$x), th.sp=0, th.env=0)  # gridding in this function there is something
    z2.sim <- ecospat.grid.clim.dyn(z2$glob, z2$glob1, data.frame(sp2.sim), R = length(z2$x), th.sp=0, th.env=0)
    
    o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = F)  # overlap between random and observed niches FAlSE because all environments
    sim.o.D <- o.i$D  # storage of overlaps
    sim.o.I <- o.i$I
    
  }
  return(c(sim.o.D, sim.o.I))
}


"niche.equivilency" <- function (z1, z2, rep, alternative = "greater", ncores = 1) 
{
  R <- length(z1$x) #100
  l <- list()
  obs.o <- ecospat.niche.overlap(z1, z2, cor = FALSE) 	#observed niche overlap, FAlSE because all environments
  if (ncores == 1) {
    sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.eq.gen, z1=z1, z2=z2)), byrow = TRUE, ncol = 2))
  }
  
  colnames(sim.o) <- c("D", "I")
  l$sim <- sim.o
  l$obs <- obs.o
  if (alternative == "greater") {
   l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) +  1)
    l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 1)
  }
  if (alternative == "lower") {
    l$p.D <- (sum(sim.o$D <= obs.o$D))/(length(sim.o$D)) ###why is it +1 
    
    l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)
  }
  return(l)
}

### For mapping densities and overlap
col_E <- rgb(1,0,0,0.2) #non-native colours
col_N <- rgb(0,0,1,0.2) #native colours

"map_overlap" <- function(y, status){
  y_E <- y[status == "E"]
  y_N <- y[status == "N"]
  dens_E <- density(y_E) 
  dens_N <- density(y_N)
  xlim <- range(dens_E$x, dens_N$x) #x axis range
  ylim <- range(0, dens_E$y, dens_N$y)
  ylim[2] <- ylim[2] + 0.1
  hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = species82[i], xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
  hist(y_N, proba = T, add = T, col = 'white', density = 10, angle = 45)
  polygon(dens_E, density = -1, col = col_E)
  polygon(dens_N, density = -1, col = col_N)
  mat <- cbind(dens_E$x, dens_E$y)
  mat <- rbind(mat, mat[1,])
  pol_E <- st_polygon(list(mat))
  mat <- cbind(dens_N$x, dens_N$y)
  mat <- rbind(mat, mat[1,])
  pol_N <- st_polygon(list(mat))
  pol_inter <- st_intersection(st_buffer(pol_E,0),  st_buffer(pol_N,0))
  plot(pol_inter, add = T, col = "grey", ylim=c(0,1))
  rug(y_N, side = 1, line = -0.15, col = col_N, tck = 0.03, lwd=2)
  rug(y_E, side = 1, line = -0.15, col = col_E, tck = 0.03, lwd=2)
}


####	####
## Calculating niche shifts for each species 
# If you do not want to run this part, you can use directly the niche_results and sample_size objects
####	####

#### Extracting occurences and calculating Overlap Niche expansion and Equivilency Test
niche_expansion_1 <- NULL
sample_size_1 <- NULL 
gridnumbers <- NULL
par(mfrow=c(1,1))
par(mar=c(5, 5, 5, 5))

#--> we can store all the density graphs in a .pdf : do not hesitate to change the parameter if this solution is not as you want
pdf("niche_plots.pdf")	# save all the overlap graphs in one pdf file
for(i in 1:length(species82)){
  ### Exctracting occurences
  species <- antmaps_occ[[species82[i]]]
  #species <- antmaps_occ[[species82[3]]]
  
  ## native data
  Native <- species[which(species$status=='N'),]
  rownames(Native) <- 1:nrow(Native)
  lat <- Native$lat 
  lon <- Native$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  N_origional <- nrow(coords)
  	#--> function below uses  nndist to calculae distance between two points
  coords <- ecospat.occ.desaggregation(xy=coords, min.dist=0.02, by=NULL) 
  N_thinned <- nrow(coords)
  index <- match(rownames(coords), rownames(Native))
  points <- SpatialPointsDataFrame(coords, data=Native[index,], proj4string = r@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(r,points, cellnumbers=T) #extracts the bioclim data for the coodrinate points we have 
  ## See how many points ecist within the same grid cell 
  grid_N <- as.data.frame(table(values[,1]))
  grid_N <- nrow( grid_N[grid_N$Var1[ grid_N$Freq > 1],])
  
  status <- points$status
  Native_df <- cbind.data.frame(coordinates(points), status, values)  #binds the bioclim data 
  Native_df <- Native_df[!duplicated(Native_df$cells),]
  Native_df <- Native_df[,-4]

  
  ## non-native data
  Exotic <- species[which(species$status=='E'),]
  rownames(Exotic) <- 1:nrow(Exotic)
  lat <- Exotic$lat
  lon <- Exotic$lon
  coords <- data.frame(x=lon,y=lat)
  E_origional <- nrow(coords)
  coords<-ecospat.occ.desaggregation(xy=coords, min.dist=0.02, by=NULL) 
  E_thinned <- nrow(coords)
  index <- match(rownames(coords), rownames(Exotic))
  points <- SpatialPointsDataFrame(coords, data=Exotic[index,], proj4string = r@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(r,points, cellnumbers=T)
  ## See how many points ecist within the same grid cell 
  grid_E <- as.data.frame(table(values[,1]))
  grid_E <- nrow(grid_E[ grid_E$Var1[ grid_E$Freq > 1],])
  spec_grid <- cbind(grid_N, grid_E, species82[i])
  colnames(spec_grid) <- c('Grid_N', 'Grid_E', "Species")
  gridnumbers <- rbind(gridnumbers, spec_grid)
  
  status <- points$status
  Exotic_df <- cbind.data.frame(coordinates(points), status, values)  #binds the bioclim data 
  Exotic_df <-  Exotic_df[!duplicated(Exotic_df$cells),]
  Exotic_df <-  Exotic_df[,-4]
  
 
  ## save number of occurences
  sizes <- cbind( N_origional, N_thinned , E_origional, E_thinned )
  rownames(sizes) <- c(species82[i])
  sample_size_1 <- rbind( sample_size_1, sizes)
  
  ## buid raw dataset for species
  #first remove all NAs
  Exotic_df <- Exotic_df[complete.cases(Exotic_df), ] 
  Native_df <- Native_df[complete.cases(Native_df), ]
  #dataset for the whole envirnoment
  Native_exotic_df <- rbind(Exotic_df, Native_df)
  status <- as.factor(Native_exotic_df$status)

  ### PCA and between class analysis
  ## perform a principal component analysis of a data frame 
  pca.env <- dudi.pca(Native_exotic_df[,4:20], scannf = F, nf = 10) ##so can first make a pca of all, look at eigen values ect. 

  ## perform between-class analysis : highlihgts the differences between groups 
  status <- as.factor(status)
  bet1 <- bca(pca.env, status, scan = FALSE, nf = 9) 
  #bet1$ratio #between-class ratio
 
  ## explore correlaitons to each bio variable to axis
  correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis. 
  correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
  rownames(correlation_to_pca) <-  c(species82[i])
  correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
  correlation_to_bio_variables <- t(correlation_to_bio_variables) 
  rownames(correlation_to_bio_variables) <-  c(species82[i])
  #ecospat.plot.contrib(bet1$co)

  ## save interesting values  
  between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
  data <-  cbind(status, between_class) #binding the two together
  #the scores for the whole environment
  scores.globclim <- data$CS1
  #scores for the presences 
  scores.sp.ex <- data$CS1[which(data$status=='E')]
  scores.sp.nat <-data$CS1[which(data$status=='N')] 
  #scores for the whole study region, same as presences as no considering whole study area 
  scores.clim.ex <- data$CS1[which(data$status=='E')]
  scores.clim.nat <-data$CS1[which(data$status=='N')] 

  ### Evaluate niche shifts
  ## Grid the whole environment 
  #correct species occupancy and whole enviornment into densities and can take climatic avaibility into account 
  grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
  #glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only, 
  grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
  
  ## Overlap
  # D overlap
  D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap. 

  # Niche expansion : delimiting niche categories and quantifying niche dynamics in analogue climates
  niche.dyn.whole <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = NA) #if intersection= NA means no intersection of analogous envirnoment
  niche.dyn.whole$dynamic.index.w 

  ## Equivilency Test
  eq.test <- NA
  eq.test <- niche.equivilency(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "lower", ncores = 1)
 
  ## save results
  row <- cbind.data.frame(D.overlap, eq.test$p.D, niche.dyn.whole$dynamic.index.w[[1]],niche.dyn.whole$dynamic.index.w[[2]], niche.dyn.whole$dynamic.index.w[[3]])
  names(row) <- c('D.overlap', 'eq.test', 'expansion', 'stability', 'Unfilling')
  rownames(row) <- c(species82[i])
  row <- cbind ( row,  correlation_to_pca,  correlation_to_bio_variables) #add the correlations of bio and pca axis to each row 
  niche_expansion_1 <- rbind(niche_expansion_1, row)

  ## plot results
  map_overlap(bet1$ls[,1], status)
  mypath3 <- file.path("overlap_graphs//",paste("density", species82[i], ".jpg", sep = "")) #path to save file 
  dev.copy(png, file=mypath3,  width = 7.7, height = 6.5, units='in', res=500)
  dev.off() #no more added to the file 
  
  
}
  dev.off()

niche_expansion_1
sample_size_1

# the data frame to be used 
	# niche_results <- niche_expansion_1[,1:5]
	# sample_size <- sample_size_1


####	####
## Evaluating dispersion groups for new species (from Bertelsmeier et al. (2017))
# If you do not want to run this part, you can use directly the groups82 object
####	####

#### Find our 82 non-native species among the 13104 species
index <- match(species82, row.names(xy_map))
species82[is.na(index)]
which(is.na(index))
index[16] <- which(row.names(xy_map) == "Monomorium latinode")
index[74] <- which(row.names(xy_map) == "Tetramorium caespitum")
xy <- xy_map[index,]

#### redo a clustering with all the species
### Add dispersion data to the original dispersion values
index <- match(row.names(xy), species241)
sum(!is.na(index))	# 25 species were not classified
new <- rbind(xy_map[species241,], xy[is.na(index),])

### Redo clustering on all species
h.group <- hclust(dist(new), "ward.D")
fac.group_1 <- as.factor(cutree(h.group, k = 4))
fac.group_1 <- as.data.frame(fac.group_1)
fac.group_1$species <- rownames(new)


# the dataframe to be used : we keep the first strategy
	# groups82 <- fac.group_1[row.names(xy),]
	# row.names(groups82) <- species82
	# groups82$species <- species82



####	####
## Calculating rao diversity in native range from Bertelsmeier et al. (2017)
# If you do not want to run this part, you can use directly the niche_results object
####	####

### Get coordinates for the middle of each bentity
sp <- coordinates(w_rworldmap) #these are the centroids of this. This is a matrix of coordinates 
sp <- SpatialPoints(sp, proj4string = w_rworldmap@proj4string) #spatial points data frame by joining the coodinates to the attributes. 
spdf <- SpatialPointsDataFrame(sp, data=w_rworldmap@data, proj4string = w_rworldmap@proj4string ) #so this is specifying to use the same system as in the worldmap one

#sp <- coordinates(map) #these are the centroids of this. This is a matrix of coordinates 
#sp <- SpatialPoints(sp, proj4string = map@proj4string) #spatial points data frame by joining the coodinates to the attributes. 
#spdf <- SpatialPointsDataFrame(sp, data=map@data, proj4string = map@proj4string ) #so this is specifying to use the same system as in the worldmap one

nr <- ne <- NULL
ns <- length(antmaps_distri)
for(i in 1:ns){
  nr <- c(nr, antmaps_distri[[i]]$gid[antmaps_distri[[i]]$category == "N"])
  ne <- c(ne, antmaps_distri[[i]]$gid[antmaps_distri[[i]]$category == "E"])
}

rownames(antmaps_spol) <- antmaps_spol$id
rownames(antmaps_spol_i) <- antmaps_spol_i$id
spol_t <- rbind(antmaps_spol, antmaps_spol_i)

lev <- as.character(spol_t$id)
nr <- factor(nr, levels = lev)
nr <- unclass(table(nr))
ne <- factor(ne, levels = lev)
ne <- unclass(table(ne))
spol_t <- cbind(spol_t, nr, ne) # sf collection
spol_t_spdf <- as(spol_t, Class = "Spatial") #wants it to be spatial polygon data frame

### Get geographic distances between each bentity
coord <- coordinates(spol_t_spdf) #map is a sparialpolygondataframe
dgeo <- matrix(0, nrow = nrow(coord), ncol = nrow(coord))
for(i in 1:nrow(coord)){
  for(j in 1:(i)){
dgeo[i,j] <- distHaversine(coord[i,], coord[j,])
  }
}
dgeo <- dgeo+t(dgeo)
dgeo <- as.data.frame(dgeo)
row.names(dgeo) <- row.names(coord)
names(dgeo) <- row.names(coord)
dgeo <- dgeo / 1000000

### Calculate Rao diversity of native range
species <- species82[1]
local <- antmaps_distri[[species]]
local <- local[which(local$category=='N'),]
df <-  as.data.frame(matrix(0, nrow=length(spol_t_spdf$gid), ncol=1))
rownames(df)  <- spol_t_spdf$gid
index <- match(local$gid, rownames(df))
df$V1[index] <- 1
colnames(df) <- species
world2 <-  df
for (i in 2:length(species82)) {
  species <- species82[i]
  local <- antmaps_distri[[species]]
  local <- local[which(local$category=='N'),]
  df <-  as.data.frame(matrix(0, nrow=length(spol_t_spdf$gid), ncol=1))
  rownames(df)  <- spol_t_spdf$gid
  index <- match(local$gid, rownames(df))
  df$V1[index] <- 1
  colnames(df) <- species
  world2<- cbind(world2, df)
  
}
colnames(world2) <- gsub(" ", "_", colnames(world2))
rao_native <- divc(world2, quasieuclid(as.dist(dgeo))) 

# add to niche_results
	# niche_results <- cbind.data.frame(niche_results, rao_native)
	# names(niche_results)[6] <- "rao_native"


####	####
## Calculating hypervolume of native niche
# This part needs long time of calcul : you can use directly the object volume 2
####	####

#### Preparing the data
### Code getting the values for all the species 
xy <- antmaps_occ[[1]][,c("lon", "lat")] #for the first species 
species <- rep(speciesnames[1], nrow(xy))
l <- lapply(antmaps_occ, nrow) # number of rows for each pccurance
index <- which(!unlist(lapply(l, is.null))) # remove nulls 
for(i in index[-1]){ #for every species appart from the first(which is already done)
  local <- antmaps_occ[[i]][,c("lon", "lat")] #get the long and lats
  xy <- rbind.data.frame(xy, local) #bind them all together 
  species <- c(species, rep(speciesnames[i], nrow(local))) #get the species names as a list
}
df <- as.data.frame(as.factor(species))
names(df) <- "species"
spdf <- SpatialPointsDataFrame(xy, df, proj = CRS(proj4string(r)))
sf <- st_as_sf(spdf) # the points and species name. 

### Format the worldclim data 
## Crop the world to the extent of occurances
e <- extent(sf) # this is the extent of ant occurances 
x1 <- -180 ; x2 <- 180 ; y1 <- round(attributes(e)$ymin)-5 ; y2 = round(attributes(e)$ymax)+5 # keep the x's the same as normal but crop the lattitudes. 
worldclim_5_crop <- crop(r, extent(x1, x2, y1, y2)) # 
 
## Downsize resolution for faster computation
worldclim_5_crop_agg <- aggregate(worldclim_5_crop, fact = 5, fun = mean)
worldclim_5_tab <- values(worldclim_5_crop)
 
# remove the nas
local <- apply(is.na(worldclim_5_tab), 1, sum)
index <- which(local == 0)
worldclim_5_tab_NA <- worldclim_5_tab[index,]

# make a pca of the whole world 
worldclim_25_pca <- dudi.pca(worldclim_5_tab_NA, nf=8, scannf=F) # pca 
#add this to the raster (so we can assign to the occurance points)
worldclim_5_tab <- as.data.frame(worldclim_5_tab)
worldclim_5_tab$pca1 <- NA
worldclim_5_tab$pca1[index] <- worldclim_25_pca$li$Axis1 # here for each of the pixels we know what the pca1 is `na or not
worldclim_5_tab$pca2 <- NA
worldclim_5_tab$pca2[index] <- worldclim_25_pca$li$Axis2 
worldclim_5_tab$pca3 <- NA
worldclim_5_tab$pca3[index] <- worldclim_25_pca$li$Axis3
worldclim_5_tab$pca4 <- NA
worldclim_5_tab$pca4[index] <- worldclim_25_pca$li$Axis4
worldclim_5_tab$pca5 <- NA
worldclim_5_tab$pca5[index] <- worldclim_25_pca$li$Axis5
worldclim_5_tab$pca6 <- NA
worldclim_5_tab$pca6[index] <- worldclim_25_pca$li$Axis6
worldclim_5_tab$pca7 <- NA
worldclim_5_tab$pca7[index] <- worldclim_25_pca$li$Axis7 
worldclim_5_tab$pca8 <- NA
worldclim_5_tab$pca8[index] <- worldclim_25_pca$li$Axis8


#### NON FIXED Bandwidth method - to calcluate optimum bandwidths
speciesnames <- as.factor(species82)
for (i in 1:length(speciesnames)) {
  ##extract the points for each species 
  species <- antmaps_occ[[as.character(speciesnames[i])]]
  ##seperate the distirbutions by native and exotic! (and then maybe by continant) 
  Native <- species[which(species$status=='N'),]
  rownames(Native) <- 1:nrow(Native)
  lat <- Native$lat 
  lon <- Native$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  N_origional <- nrow(coords)
  ##function below uses  nndist to calculae distance ebtween two points
  coords<-ecospat.occ.desaggregation(xy=coords, min.dist=0.02, by=NULL) ###assume same as worldclim  (arc rotations are the unit) so 0.5 is 1km
  N_thinned <- nrow(coords)
  index <- match(rownames(coords), rownames(Native))
  points <- SpatialPointsDataFrame(coords, data=Native[index,], proj4string = r@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(worldclim_5_crop, points, cellnumbers=T) #extracts the bioclim data for the coodrinate points we have 
  Native_df <- cbind.data.frame(coordinates(points),  values)  #binds the bioclim data 
  #plot(r[[1]]) #just to check if the points plot as expected
  
  #matches pca values to the occurance points 
  Native_df <- Native_df[complete.cases(Native_df), ]
  index <- match(as.character(Native_df$cells), rownames(worldclim_5_tab))
  Native_df$pca1 <- worldclim_5_tab$pca1[index]
  Native_df$pca2 <- worldclim_5_tab$pca2[index]
  Native_df$pca3 <- worldclim_5_tab$pca3[index]
  Native_df$pca4 <- worldclim_5_tab$pca4[index]
  Native_df$pca5 <- worldclim_5_tab$pca5[index]
  Native_df$pca6 <- worldclim_5_tab$pca6[index]
  Native_df$pca7 <- worldclim_5_tab$pca7[index]
  Native_df$pca8 <- worldclim_5_tab$pca8[index]
  
  #make a hypervolume, as based on the same pca axis then it should be comparable between species 
  #estimate bandwidth for each species, then use the highest one. 
  hyper_world <- hypervolume_gaussian(Native_df[,21:25], chunk.size=500) #gaussian method is better as the box method is flat therefore has errors 
  
  #dataset with badwidth details 
  bandwidth_sp <- hyper_world@Parameters$kde.bandwidth # this gives us the bandwidth. 
  bandwidth_sp <- as.data.frame(bandwidth_sp)
  bandwidth_sp <- t(bandwidth_sp)
  rownames(bandwidth_sp) <- as.character(speciesnames[i])
  bandwidth_details <- rbind(bandwidth_details, bandwidth_sp)
  
  #dataset with the volume details 
  volume_sp <- get_volume(hyper_world)
  volume_sp <- as.data.frame(volume_sp)
  rownames(volume_sp) <- as.character(speciesnames[i])
  volume <- rbind(volume, volume_sp)
}
######

#remove the outliers ! 
library(outliers)
bandwidth_details <- as.data.frame(bandwidth_details)
while (grubbs.test(bandwidth_details$pca1)$p.value <0.05) {
  bandwidth_details$pca1[which(bandwidth_details$pca1==outlier(bandwidth_details$pca1))] <- NA
}
max(na.omit(bandwidth_details$pca1)) # 1.56

while (grubbs.test(bandwidth_details$pca2)$p.value <0.05) {
  bandwidth_details$pca2[which(bandwidth_details$pca2==outlier(bandwidth_details$pca2))] <- NA
}
max(na.omit(bandwidth_details$pca2)) #  2.19

while (grubbs.test(bandwidth_details$pca3)$p.value <0.05) {
  bandwidth_details$pca3[which(bandwidth_details$pca3==outlier(bandwidth_details$pca3))] <- NA
}
max(na.omit(bandwidth_details$pca3)) # 1.5

while (grubbs.test(bandwidth_details$pca4)$p.value <0.05) {
  bandwidth_details$pca4[which(bandwidth_details$pca4==outlier(bandwidth_details$pca4))] <- NA
}
max(na.omit(bandwidth_details$pca4)) # 0.67

while (grubbs.test(bandwidth_details$pca5)$p.value <0.05) {
  bandwidth_details$pca5[which(bandwidth_details$pca5==outlier(bandwidth_details$pca5))] <- NA
}
max(na.omit(bandwidth_details$pca5)) #  1.05

### ##
##### #
### FIXED BANDWIDTH METHOD 
#### #
### ##
bwdths <- c(1.56, 2.19, 1.50, 0.67, 1.05)

## using fixed bandwidth of =previous section
#LOOPS Native
bandwidth_details2 <- NULL 
volume2_1 <- NULL 
for (i in 1:length(speciesnames)) {
  ## extract the points for each species 
  species <- antmaps_occ[[as.character(speciesnames[i])]]
  ## seperate the distirbutions by native and exotic! (and then maybe by continant) 
  Native <- species[which(species$status=='N'),]
  rownames(Native) <- 1:nrow(Native)
  lat <- Native$lat 
  lon <- Native$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  N_origional <- nrow(coords)
  ## function below uses  nndist to calculae distance ebtween two points
  coords<-ecospat.occ.desaggregation(xy=coords, min.dist=0.02, by=NULL) ###assume same as worldclim  (arc rotations are the unit) so 0.5 is 1km
  N_thinned <- nrow(coords)
  index <- match(rownames(coords), rownames(Native))
  points <- SpatialPointsDataFrame(coords, data=Native[index,], proj4string = r@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(worldclim_5_crop, points, cellnumbers=T) #extracts the bioclim data for the coodrinate points we have 
  Native_df <- cbind.data.frame(coordinates(points),  values)  #binds the bioclim data 
  
  # matches pca values to the occurance points 
  index <- match(as.character(Native_df$cells), rownames(worldclim_5_tab))
  Native_df$pca1 <- worldclim_5_tab$pca1[index]
  Native_df$pca2 <- worldclim_5_tab$pca2[index]
  Native_df$pca3 <- worldclim_5_tab$pca3[index]
  Native_df$pca4 <- worldclim_5_tab$pca4[index]
  Native_df$pca5 <- worldclim_5_tab$pca5[index]
  Native_df$pca6 <- worldclim_5_tab$pca6[index]
  Native_df$pca7 <- worldclim_5_tab$pca7[index]
  Native_df$pca8 <- worldclim_5_tab$pca8[index]
  Native_df <- Native_df[complete.cases(Native_df[,21:25]), ]
  # make a hypervolume, as based on the same pca axis then it should be comparable between species 
  # estimate bandwidth for each species, then use the highest one. 
  hyper_world <- hypervolume_gaussian(Native_df[,21:25], chunk.size=500, kde.bandwidth=bwdths) #gaussian method is better as the box method is flat therefore has errors 
  
  # dataset with badwidth details 
  bandwidth_sp <- hyper_world@Parameters$kde.bandwidth # this gives us the bandwidth. 
  bandwidth_sp <- as.data.frame(bandwidth_sp)
  bandwidth_sp <- t(bandwidth_sp)
  rownames(bandwidth_sp) <- as.character(speciesnames[i])
  bandwidth_details2 <- rbind(bandwidth_details2, bandwidth_sp)
  
  # dataset with the volume details 
  volume_sp <- get_volume(hyper_world)
  volume_sp <- as.data.frame(volume_sp)
  rownames(volume_sp) <- as.character(speciesnames[i])
  volume2_1 <- rbind(volume2_1, volume_sp)
}

# add to niche_results
	# niche_results <- cbind.data.frame(niche_results, volume2_1)
	# names(niche_results)[7] <- "hypervolume_native"


####	####
## Statistics
####	####

#### Preparing data
niche_results$X <- gsub(" ", "_", row.names(niche_results))

## Group 0-10 and 10-100% expansion
niche_results$expansion_sig  <- ifelse(niche_results$expansion > 0.1, 'expansion', niche_results$expansion)
niche_results$expansion_sig  <- ifelse (niche_results$expansion_sig  <= 0.1, 'noexpansion', niche_results$expansion_sig)

## Overview of data
par(mar = c(4, 4, 4, 4))
hist(niche_results$D.overlap, col='grey', breaks=20, xlab='D Overlap', main='D Overlap of Native and Invasive Ranges') 
hist(niche_results$expansion, col='grey', breaks=20, xlab='Percentage of Exotic Range Expansion', main='Expansions Observed in Invasive Range') 
shapiro.test(niche_results$D.overlap) #significant - need to use non-paremetric tests
shapiro.test(niche_results$expansion) #significant - need to use non-paremetric tests
# mean and standard deviation
sumstat <- niche_results %>%
# Select and rename five variables 
dplyr::select(`D Overlap` = D.overlap,`Expansion (%)` = expansion) %>%
# Find the mean, st. dev., min, and max for each variable 
summarise_each(funs(mean, sd, min, max, median)) %>%  # mean, st. dev., min, and max for each variable 
# Move summary stats to columns
gather(key, value, everything()) %>% 
separate(key, into = c("variable", "stat"), sep = "_") %>%
spread(stat, value) %>%
# Set order of summary statistics 
dplyr::select(variable, median, mean, sd, min, max) %>%
# Round all numeric variables to one decimal point
mutate_each(funs(round(., 2)), -variable)
sumstat

## Sample size mean and median
sample_size <- as.data.frame(sample_size)
sampletotal<- sample_size$N_origional + sample_size$E_origional
#mean(sampletotal)
#std.error(sampletotal)
median(sampletotal)
min(sampletotal)
max(sampletotal)

#Equivilency data. As there are multiple species we are comparing between, need to adjust p-values
niche_results$eq.testBH  <-    p.adjust(niche_results$eq.test ,     method = "BH")


####	####
## Figure 1
####	####


####Then the map 
#For each of the four example species 
spe <-"Tetramorium bicarinatum" #A
#spe <- "Strumigenys margaritae" #B
#spe <-"Amblyopone australis" #C
#spe <- "Monomorium pharaonis" #D
### Occurences data
xy <- antmaps_occ[[spe]][,c("lon", "lat")]
colnames(xy) <- c('x', 'y')
xy<-ecospat.occ.desaggregation(xy=xy, min.dist=0.03, by=NULL) ###assume same as worldclim  (arc rotations are the unit) so 0.5 is 1km
species2 <- rep(spe, nrow(xy))
index <- match(rownames(xy), rownames( antmaps_occ[[spe]][,c("lon", "lat")]))
status2 <- antmaps_occ[[spe]]$status[index]
l <- lapply(antmaps_occ, nrow)
index2 <- which(!unlist(lapply(l, is.null)))

df <- cbind.data.frame(as.factor(species2), as.factor(status2))
names(df) <- c("species", "status")

spdf <- SpatialPointsDataFrame(xy, df, proj = r@crs)
sf <- st_as_sf(spdf)
crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
sf <- lwgeom::st_transform_proj(sf, crs = crs_wintri)


#plotting the native and exotic ranges
sf_N <- sf[sf$status == "N",] # native range
sf_E <- sf[sf$status == "E",] #exotic range 

species <- antmaps_occ[[spe]]
Native <- species[which(species$status=='N'|species$status=='E'),]
Native <- species[which(species$status %in% c('N','E')),]
lat <- Native$lat 
lon <- Native$lon
coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 

#creating the map
w_rworldmap_crop <- raster::crop(w_rworldmap, extent(-170, 170, -65, 90))
w_rworldmap_sf <- st_as_sf(w_rworldmap_crop)
w_rworldmap_sf_wintri <- lwgeom::st_transform_proj(w_rworldmap_sf, crs='+proj=wintri')

#plotting
par(mar=c(0,0,0,0))
plot <- plot(st_geometry(w_rworldmap_sf_wintri),  ylim=c( -65, 90),  col='grey80', border='grey40', lwd=0.6)
plot <- plot + plot(sf_N$geom, add = T, pch = "+", col='royalblue1', cex=1.7)
plot <- plot + plot(sf_E$geom, add = T, pch = "+", col = "red", cex=1.7)
dev.off()
#1500, 1150






####	####
## Figure 2
####	####

#### Alien Vs Invasive
niche_results$type <- status82      

### Effect on D overlap 
kruskal.test(D.overlap ~ type, data = niche_results) # SIGNIFICANT
## plot it (fig 2.A)
boxA <- ggplot( niche_results , aes(x=type, y=D.overlap, fill=type)) +
 scale_fill_manual(values=c("#c3dfe0",  '#2b7570'), name='Expansion Groups', labels=c('Alien', 'Invasive')) +
 geom_boxplot(color='black', outlier.shape = NA) + 
 geom_jitter(shape=19, position=position_jitter(0.2), size=0.7) + 
 xlab('') +  ylab('D overlap') + 
 theme_classic()  + 
 scale_x_discrete(breaks=c("A","I"),  labels=c("Alien", "Invasive")) + 
 theme(legend.position = 'none') + 
 theme(text = element_text(family = "Arial", size=10)) + 
 theme(axis.text.x=element_text(colour="black")) +
 annotate("text", x = 1.5, y =0.9 , label = "*", size=4) +  
 geom_segment(aes(x = 1, y = 0.88, xend = 2, yend = 0.88,), colour='black')
boxA

### Effect on expansion
## Need to create the Dataset 
expansion_categories <- as.data.frame(cbind(niche_results$expansion_sig, niche_results$type))
expansion_categories <- expansion_categories[complete.cases(expansion_categories),]
expansion_categories$V1 <- as.factor(expansion_categories$V1)
expansion_categories$V2 <- as.factor(expansion_categories$V2)
levels(expansion_categories$V2) <- c("A", "I")
# percentage of alien species spread
expansion_categories_alien <- expansion_categories[which(expansion_categories$V2=='A'),]
expansion_test_data_alien  <- expansion_categories_alien  %>% 
                               dplyr::count(V1) %>%  
                               mutate(perc = n / nrow( x =expansion_categories_alien))
expansion_test_data_alien$V2 <- rep('A', nrow(expansion_test_data_alien))
expansion_test_data_alien <- expansion_test_data_alien %>% arrange(desc(V1)) 
# percentage of invasive species spread
expansion_categories_inv <- expansion_categories[which(expansion_categories$V2=='I'),]
expansion_categories_inv  %>%
 dplyr::count(V1 )%>%  
 dplyr::mutate(perc = n / nrow( expansion_categories_inv)) -> expansion_test_data_inv
expansion_test_data_inv$V2 <- rep('I', nrow(expansion_test_data_inv))
expansion_test_data_inv <- expansion_test_data_inv %>% arrange(desc(V1)) 
# combine them to plot them both on a chart  
expansion_test_data <- rbind(expansion_test_data_inv, expansion_test_data_alien)
expansion_test_data <- as.data.frame(expansion_test_data)

## Plot it (fig 2.B)
colours <- c( 'coral3', 'gray83')
supp.labs <- c('Alien', 'Invasive')
names(supp.labs) <- c("A", "I")
alienpieA <- ggplot(expansion_test_data, aes(x = "", y = perc, fill = V1)) +
 geom_bar(width = 1, stat = "identity", color = "white") +
 coord_polar("y", start = 0)+ facet_wrap(~ V2, labeller=labeller(V2=supp.labs)) +
 geom_text(aes(label = percent(round(perc*100, digits=1)/100)), position = position_fill(vjust = 0.5), size=3.5) +
 scale_fill_manual(values = colours)  + 
 theme_void() + theme(plot.title = element_text(hjust = 0.5))  + 
 labs(fill = "expansion") + 
 theme(text = element_text(family = "Arial", size=10))
alienpieA 

## Statistical test  
expansion_test <- expansion_test_data[,c(1, 2, 4)] %>%
 dcast(V2~V1, value.var='n') %>%
 as.data.frame()
expansion_test <- expansion_test[,c(2, 3)]
chisq.test(expansion_test) ##not significant


#### Global Geographic Expansion groups (second test of invasivness)  
### Preparing data
names(groups82)[1] <- "fac.group"
groups82$species <- as.factor(groups82$species)
groups82$species  <- gsub(" ", "_", groups82$species )
         
index <- match(niche_results$X, groups82$species) 
niche_results$categories <- groups82$fac.group[index]
         
# we dont want 1's
index <- which(niche_results$categories==1)
niche_results$categories[index] <- 2
groups82_comparison <- niche_results[complete.cases(niche_results$categories), ]
 
#

regional  <- niche_results[which(niche_results$categories==2),]
length(which(regional$type=='I'))
Transcont  <- niche_results[which(niche_results$categories==3),]
length(which(Transcont$type=='I'))
glob <- niche_results[which(niche_results$categories==4),]
length(which(glob$type=='I'))
       
### Effect on D overlap      
# Stats
kruskal.test(D.overlap ~ as.factor(categories), data =  groups82_comparison) # this is less than 0.05
pairwise.wilcox.test(groups82_comparison$D.overlap, as.factor(groups82_comparison$categories), p.adjust.method="BH") 

## Plot it (fig2.C)
boxB <- ggplot(groups82_comparison , aes(x=as.factor(categories), y=D.overlap, fill=as.factor(categories))) +
 geom_boxplot(color='black', outlier.shape = NA) + 
 geom_jitter(shape=19, position=position_jitter(0.2), size=0.7) + 
 scale_fill_manual(values=c("#A8DADC", '#028090', '#05668d'),name='Dispersion Groups', labels=c('Regional', 'Transcontinental', 'Global')) +
 xlab('') + ylab('D overlap') + 
 theme_classic() + theme(legend.position='none') + 
 scale_x_discrete(breaks=c("2","3", "4"),  
 labels=c('Regional', 'Transcontinental', 'Global')) +
 theme(text = element_text(family = "Arial", size=10)) +
 theme(axis.text.x=element_text(colour="black")) +
 annotate("text", x = 2, y =0.91 , label = "***", size=4) +  
 geom_segment(aes(x = 1, y = 0.9, xend = 3, yend = 0.9), colour='black') +
 annotate("text", x = 1.5, y =0.71 , label = "***", size=4) + 
 geom_segment(aes(x = 1, y = 0.7, xend = 2, yend = 0.7), colour="black") + 
 annotate("text", x = 2.5, y =0.84 , label = "*", size=4) +  
 geom_segment(aes(x = 2, y = 0.82, xend = 3, yend = 0.82), colour="black")
boxB

###  Effect on expansion
## Create Dataset of percentages of each category
expansion_categories <- as.data.frame(cbind(groups82_comparison$expansion_sig, groups82_comparison$categories ))
expansion_categories <- expansion_categories[complete.cases(expansion_categories),]
expansion_categories$V1 <- as.factor(expansion_categories$V1)
expansion_categories$V2 <- as.factor(expansion_categories$V2)
#percentage of alien species spread
expansion_test_data_2 <- expansion_categories[which(expansion_categories$V2=='2'),]  %>%
 dplyr::count(V1) %>%  
 mutate(perc = n / nrow(expansion_categories[which(expansion_categories$V2=='2'),]))#percentages. 
expansion_test_data_2$V2 <- rep('2', nrow(expansion_test_data_2))
expansion_test_data_2 <- expansion_test_data_2 %>%
 arrange(desc(V1)) 
#percentage of invasive species spread
expansion_test_data_3 <- expansion_categories[which(expansion_categories$V2=='3'),]  %>%
 dplyr::count(V1 )%>%  
 dplyr::mutate(perc = n / nrow(expansion_categories[which(expansion_categories$V2=='3'),]))  
expansion_test_data_3$V2 <- rep('3', nrow(expansion_test_data_3))
expansion_test_data_3 <- expansion_test_data_3 %>%
 arrange(desc(V1)) 
#combine them to plot them both on a chart  
expansion_test_data_4 <- expansion_categories[which(expansion_categories$V2=='4'),] %>%
 dplyr::count(V1 )%>%  
 dplyr::mutate(perc = n / nrow( expansion_categories[which(expansion_categories$V2=='4'),])) 
expansion_test_data_4$V2 <- rep('4', nrow(expansion_test_data_4))
expansion_test_data_4 <- expansion_test_data_4 %>%
 arrange(desc(V1)) 
# combine
 expansion_test_data <- rbind(expansion_test_data_2, expansion_test_data_3, expansion_test_data_4)
         
## Plot it
colours <- c( 'coral3', 'grey')
supp.labs <- c('Regional', 'Transcontinental', 'Global')
names(supp.labs) <- c("2", "3", "4")
pieA.1 <- ggplot(expansion_test_data, aes(x = "", y = perc, fill = V1)) +
 geom_bar(width = 1, stat = "identity", color = "white") +
 coord_polar("y", start = 0)+ facet_wrap(~ V2, labeller=labeller(V2=supp.labs)) +
 geom_text(aes(label = scales::percent(round(perc*100, digits=1)/100)), position = position_fill(vjust = 0.5), size=3.5) +
 scale_fill_manual(values = colours)  + 
 theme_void() +   theme(plot.title = element_text(hjust = 1)) + 
 theme(text = element_text(family = "Arial", size=10)) +
 theme(legend.position = 'none') 
pieA.1

# Statistical test
expansion_test <- expansion_test_data[,c(1, 2, 4)]
expansion_test <- dcast(expansion_test,V2~V1, value.var='n')
expansion_test <- as.data.frame(expansion_test)
expansion_test <- expansion_test[,c(2, 3)]
chisq.test(expansion_test)  ###Significant


#### Figure 2
Fig2 <- ggarrange(boxA, alienpieA, boxB, pieA.1, 
                  labels = c("a", "b", 'c', 'd'),  font.label = list(size = 10, font='Arial'),
                  ncol = 2, nrow = 2, widths=c(1, 1.25))
ggsave(filename = "fig2.pdf", plot=Fig2, width=7.2, height= 5, units = "in", device='pdf')



####	####
## Figure 3
####	####

#### Native range characteristics 
### Expansion
# rao_native$X <- rownames(rao_native)
# rao_native$X <- gsub(" ", "_", rao_native$X)
# index <- match(niche_results$X, rao_native$X)
niche_results$raodiversity_native <- niche_results$rao_native
niche_results$expansion_sig <- factor(niche_results$expansion_sig, levels=c('noexpansion', 'expansion'))
 
# Statistical test
kruskal.test(niche_results$raodiversity_native, niche_results$expansion_sig)

# Plot
box_rao <-  ggplot(niche_results , aes(x=expansion_sig, y=raodiversity_native, fill=expansion_sig))+
 scale_fill_manual(values=c('grey', '#8E3951')) + 
 geom_boxplot(color='black', outlier.shape = NA) + 
 geom_jitter(shape=16, position=position_jitter(0.1), size=0.7) + 
 scale_x_discrete(breaks=c("noexpansion","expansion"), labels=c('0-10% ', '10-100%')) +
 ylab('Native spatial \nRao diversity') + 
 xlab('Expansion (%)') + theme_classic() + theme(legend.position='none') + 
 theme(text = element_text(family = "Arial", size=10)) + 
 theme(axis.text.x=element_text(colour="black")) +
 annotate("text", x = 1.5, y =58 , label = "***", size=4) +  
 geom_segment(aes(x = 1, y = 57, xend = 2, yend = 57), colour="black") 
box_rao   
         
### D overlap
# Statistical test
cor.test(niche_results$D.overlap, niche_results$raodiversity_native, method = 'spearman') #native range

# Plot
scat_rao <- ggscatter(niche_results, x = 'raodiversity_native', y = 'D.overlap', 
 add = "reg.line", conf.int = TRUE, 
 cor.coef = TRUE, cor.method = "kendall", size=0.6, cor.coef.size=3.3, font.label=c(10, 'plain', 'black'), font.family='Arial') +
 theme(text = element_text(family = "Arial", size=10)) +
 xlab('Native spatial Rao diversity') +  ylab(expression("D Overlap"))  + 
 labs(color = "Percentage of \nnative niche")
scat_rao     
     
                
#### Hypervolume 
### Expansion
# volume2$X <- gsub(" ", "_", volume2$X)
# index <- match(niche_results$X, volume2$X)
niche_results$volume_native_fixed <- niche_results$hypervolume_native
         
# Staistical test                    
kruskal.test( niche_results$volume_native_fixed, niche_results$expansion_sig)
        
# Plot
box_hyp <- ggplot(niche_results , aes(x=expansion_sig, y=volume_native_fixed, fill=expansion_sig))+
 scale_fill_manual(values=c('grey',  '#8E3951')) + 
 geom_boxplot(color='black', outlier.shape = NA) + 
 geom_jitter(shape=16, position=position_jitter(0.1), size=0.7) + 
 scale_x_discrete(breaks=c("noexpansion","expansion"),  labels=c('0-10%', '10-100%')) +
 ylab('Native hypervolume') + 
 xlab('Expansion (%)') + theme_classic() + theme(legend.position='none') + 
 theme(text = element_text(family = "Arial", size=10)) + 
 theme(axis.text.x=element_text(colour="black")) +
 annotate("text", x = 1.5, y =38000 , label = "***", size=4) +  
 geom_segment(aes(x = 1, y = 37100, xend = 2, yend = 37100), colour="black") 
box_hyp
         
### D overlap
# Statistical test  
cor.test(niche_results$volume_native_fixed, niche_results$D.overlap , method = 'kendall') #significant

# Plot
scat_hyp <-ggscatter(niche_results, x = 'volume_native_fixed', y = 'D.overlap', 
 add = "reg.line", conf.int = TRUE, 
 cor.coef = TRUE, cor.method = "kendall", size=0.6, cor.coef.size=3.3, font.label=c(10, 'plain', 'black'), font.family='Arial') +
 scale_x_continuous(breaks = pretty(niche_results$volume_native_fixed, n = 3)) +
 theme(text = element_text(family = "Arial", size=10)) +
 xlab('Native hypervolume') +  ylab(expression("D Overlap"))  + 
 labs(color = "Percentage of \nnative niche")
scat_hyp       
         
### Check if there is an effect of invasiveness on native characteristics 
## IUCN invasivness
kruskal.test(niche_results$raodiversity_native,  niche_results$type)
kruskal.test(niche_results$volume_native_fixed,  niche_results$type)
## Global spread 
kruskal.test(niche_results$volume_native_fixed,  as.factor(niche_results$categories))

kruskal.test(niche_results$raodiversity_native,  as.factor(niche_results$categories))

dunn.test(niche_results$raodiversity_native,  as.factor(niche_results$categories), method='bonferroni')

pairwise.wilcox.test(niche_results$raodiversity_native, as.factor(niche_results$categories), p.adjust.method="bonf") 
pairwise.wilcox.test(niche_results$raodiversity_native, as.factor(niche_results$categories), p.adjust.method="BH") 

           
#### Figure 3
Fig3 <- ggarrange(box_rao, scat_rao, box_hyp, scat_hyp, 
                   labels = c("a", "b", "c", "d"),  font.label = list(size = 10, font='Arial'),
                   ncol = 2, nrow = 2)
#ggsave(filename = "plot3_DE.png", Fig3, width=4.5, height= 4, dpi = 300, units = "in", device='png')
ggsave(filename = "Fig3.pdf", plot=Fig3, width=7, height= 5, units = "in", device='pdf')
dev.off()
