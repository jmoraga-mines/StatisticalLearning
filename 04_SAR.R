if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

# pacman::p_load(sf, tidyverse, terra)
require(sf)
require(terra)
require(tidyverse)
require(dplyr)
require(raster)

# Create circular weight matrix
circular_weight <- function(radius, crs=32611, fillNA=FALSE, cell_value=NA, 
                            zero_center=FALSE) {
  empty_raster <- raster::raster(ncols=21, nrows=21, xmn=0, res=c(1, 1), crs=crs)
  Weight_matrix <- raster::focalWeight(empty_raster, radius, "circle", fillNA = fillNA)
  if(!is.na(cell_value)){
    center_value <- as.matrix(Weight_matrix)[radius+1, radius+1]
    Weight_matrix <- Weight_matrix*cell_value/center_value
  }
  if(zero_center){
    Weight_matrix[radius+1, radius+1] = 0
  }
  return(Weight_matrix)
}



### Constants


dir1 = "e:/geostats"
dir2 = "/store02/geostats"
if(dir.exists(dir1)) {
  dir_base = dir1
  } else if(dir.exists(dir2)){
  dir_base = dir2
  }
new_brady_file <- file.path(dir_base, "new_brady_som")
new_brady_som <- raster::stack(new_brady_file)
brady_df <- raster::as.data.frame(new_brady_som, xy = T)
rm(new_brady_som)

# pacman::p_load(spdep)
require(spdep)
require(spatialreg)

# distance_nb <- dnearneigh(brady_df, 0, 30)

###############################
require(caret)
training.samples <- if_else(is.na(brady_df$Geothermal), 
                            0, brady_df$Geothermal, 0) %>%
  createDataPartition(p=0.1, list=FALSE)
brady_df <- brady_df[training.samples, ]
coordinates(brady_df) <- ~x+y
brady_df <- st_as_sf(brady_df)

attach(brady_df)
###################################################################

#load(file.path(dir_base, "distance_listw_3-18.RData"))
Sys.time(); system.time(distance_nb_3 <- (knearneigh(brady_df, k=4))); Sys.time()
distance_nb_3 <- knn2nb(distance_nb_3 )
distance_listw_3 <- nb2listw(distance_nb_3)
sar.chi<-lagsarlm(Geothermal ~ Temperature+Faults+Chalcedony+Hematite+Epsomite, 
                  data=brady_df, distance_listw_3, zero.policy = T)


###################################################################
Sys.time(); system.time(distance_nb_3 <- dnearneigh(brady_df, 0, 3)); Sys.time()
Sys.time(); system.time(distance_nb_6 <- dnearneigh(brady_df, 0, 6)); Sys.time()
Sys.time(); system.time(distance_nb_9 <- dnearneigh(brady_df, 0, 9)); Sys.time()
Sys.time(); system.time(distance_nb_12 <- dnearneigh(brady_df, 0, 12)); Sys.time()
Sys.time(); system.time(distance_nb_15 <- dnearneigh(brady_df, 0, 15)); Sys.time()
Sys.time(); system.time(distance_nb_18 <- dnearneigh(brady_df, 0, 18)); Sys.time()
save(list=grep("distance_nb", ls(), value=T), file="e:/geostats/distance_nb.RData")
Sys.time(); system.time(distance_listw_3 <- nb2listw(distance_nb_3)); Sys.time()
Sys.time(); system.time(distance_listw_6 <- nb2listw(distance_nb_6)); Sys.time()
Sys.time(); system.time(distance_listw_9 <- nb2listw(distance_nb_9)); Sys.time()
Sys.time(); system.time(distance_listw_12 <- nb2listw(distance_nb_12)); Sys.time()
Sys.time(); system.time(distance_listw_15 <- nb2listw(distance_nb_15)); Sys.time()
Sys.time(); system.time(distance_listw_18 <- nb2listw(distance_nb_18)); Sys.time()
save(list=grep("distance_listw", ls(), value=T), 
     file=file.path(dir_base, "distance_list2_02.RData"))

moran_Temp_3 <- moran.test(Temperature, distance_listw_3, na.action = na.omit)
moran_Temp_6 <- moran.test(Temperature, distance_listw_6, na.action = na.omit)
moran_Temp_9 <- moran.test(Temperature, distance_listw_9, na.action = na.omit)
moran_Temp_12 <- moran.test(Temperature, distance_listw_12, na.action = na.omit)
moran_Temp_15 <- moran.test(Temperature, distance_listw_15, na.action = na.omit)
moran_Temp_18 <- moran.test(Temperature, distance_listw_18, na.action = na.omit)

moran_Chalcedony_3 <- moran.test(Chalcedony, distance_listw_3, na.action = na.omit)
moran_Chalcedony_6 <- moran.test(Chalcedony, distance_listw_6, na.action = na.omit)
moran_Chalcedony_9 <- moran.test(Chalcedony, distance_listw_9, na.action = na.omit)
moran_Chalcedony_12 <- moran.test(Chalcedony, distance_listw_12, na.action = na.omit)
moran_Chalcedony_15 <- moran.test(Chalcedony, distance_listw_15, na.action = na.omit)
moran_Chalcedony_18 <- moran.test(Chalcedony, distance_listw_18, na.action = na.omit)


moran_Chalcedony_3$statistic
moran_Chalcedony_6$statistic
moran_Chalcedony_9$statistic
moran_Chalcedony_12$statistic
moran_Chalcedony_15$statistic
moran_Chalcedony_18$statistic
moran_Chalcedony_21$statistic
moran_Chalcedony_24$statistic
moran_Chalcedony_27$statistic
moran_Chalcedony_30$statistic

#### Lag 21
Sys.time(); system.time(distance_nb_21 <- dnearneigh(brady_df, 0, 21)); Sys.time()
save(distance_nb_21, file=file.path(dir_base, "distance_nb_21.RData"))
Sys.time(); system.time(distance_listw_21 <- nb2listw(distance_nb_21)); Sys.time()
if (exists("distance_listw_21")) rm(distance_nb_21)
save(distance_listw_21, file=file.path(dir_base, "distance_listw_21.RData"))
moran_Chalcedony_21 <- moran.test(Chalcedony, distance_listw_21, na.action = na.omit)

#### Lag 24
Sys.time(); system.time(distance_nb_24 <- dnearneigh(brady_df, 0, 24)); Sys.time()
save(distance_nb_24, file=file.path(dir_base, "distance_nb_24.RData"))
Sys.time(); system.time(distance_listw_24 <- nb2listw(distance_nb_24)); Sys.time()
if (exists("distance_listw_24")) rm(distance_nb_24)
save(distance_listw_24, file=file.path(dir_base, "distance_listw_24.RData"))
moran_Chalcedony_24 <- moran.test(Chalcedony, distance_listw_24, na.action = na.omit)

#### Lag 27
Sys.time(); system.time(distance_nb_27 <- dnearneigh(brady_df, 0, 27)); Sys.time()
save(distance_nb_24, file=file.path(dir_base, "distance_nb_27.RData"))
Sys.time(); system.time(distance_listw_27 <- nb2listw(distance_nb_27)); Sys.time()
if (exists("distance_listw_27")) rm(distance_nb_27)
save(distance_listw_27, file=file.path(dir_base, "distance_listw_27.RData"))
moran_Chalcedony_27 <- moran.test(Chalcedony, 
                                  distance_listw_27, 
                                  na.action = na.omit)

#### Lag 30
Sys.time(); system.time(distance_nb_30 <- dnearneigh(brady_df, 0, 30)); Sys.time()
save(distance_nb_30, 
     file=file.path(dir_base, "distance_nb_30.RData"))
Sys.time(); system.time(distance_listw_30 <- nb2listw(distance_nb_30)); Sys.time()
if (exists("distance_listw_30")) rm(distance_nb_30)
save(distance_listw_30, 
     file=file.path(dir_base, "distance_listw_30.RData"))
moran_Chalcedony_30 <- moran.test(Chalcedony, 
                                  distance_listw_30, 
                                  na.action = na.omit)
##########################
brady_ml_10 = MoranLocal(new_brady_som[["Geothermal"]], 
                         w=circular_weight(10, zero_center = T, cell_value = 1))

for (i in seq(from=10,to=100,by=10)){
  cat(paste0("Lag ", i, "\n"))
  print(Moran(new_brady_som[["Chalcedony"]], w=circular_weight(i, zero_center = T, cell_value = 1)))
}


Moran(new_brady_som[["Geothermal"]], w=circular_weight(50, zero_center = T, cell_value = 1))


cat("Faults\n")
for (i in seq(from=1,to=20,by=1)){
  cat(paste0("Lag ", i, ": "))
  cat(Moran(new_brady_som[["Faults"]], w=circular_weight(i, zero_center = T, cell_value = 1)))
  cat("\n")
}
