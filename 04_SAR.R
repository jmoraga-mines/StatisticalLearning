if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, tidyverse, terra)

new_brady_som <- raster::stack("e:/geostats/new_brady_som")
brady_df <- raster::as.data.frame(new_brady_som, xy = T)

pacman::p_load(spdep)
coordinates(brady_df) <- ~x+y
brady_df <- st_as_sf(brady_df)
distance_nb <- dnearneigh(brady_df, 0, 30)

sys.time(); system.time(distance_nb_3 <- dnearneigh(brady_df, 0, 3)); Sys.time()
sys.time(); system.time(distance_nb_6 <- dnearneigh(brady_df, 0, 6)); Sys.time()
sys.time(); system.time(distance_nb_9 <- dnearneigh(brady_df, 0, 9)); Sys.time()
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
save(list=grep("distance_listw", ls(), value=T), file="e:/geostats/distance_nb.RData")
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

Sys.time(); system.time(distance_nb_21 <- dnearneigh(brady_df, 0, 21)); Sys.time()
Sys.time(); system.time(distance_listw_21 <- nb2listw(distance_nb_21)); Sys.time()
Sys.time(); system.time(distance_nb_24 <- dnearneigh(brady_df, 0, 24)); Sys.time()
Sys.time(); system.time(distance_listw_24 <- nb2listw(distance_nb_24)); Sys.time()
Sys.time(); system.time(distance_nb_27 <- dnearneigh(brady_df, 0, 27)); Sys.time()
Sys.time(); system.time(distance_listw_27 <- nb2listw(distance_nb_27)); Sys.time()
Sys.time(); system.time(distance_nb_30 <- dnearneigh(brady_df, 0, 30)); Sys.time()
Sys.time(); system.time(distance_listw_30 <- nb2listw(distance_nb_30)); Sys.time()
moran_Chalcedony_21 <- moran.test(Chalcedony, distance_listw_21, na.action = na.omit)
moran_Chalcedony_24 <- moran.test(Chalcedony, distance_listw_24, na.action = na.omit)
moran_Chalcedony_27 <- moran.test(Chalcedony, distance_listw_27, na.action = na.omit)
moran_Chalcedony_30 <- moran.test(Chalcedony, distance_listw_30, na.action = na.omit)
