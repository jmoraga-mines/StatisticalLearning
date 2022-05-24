if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, tmap, dplyr)

# Reads 
# Original files in: /store03/geodata/RawData/Mahmut/SAR/Results/PSI/Brady_and_Desert/2017-2019/72Images
SarprozFile <- "d:/Geocomputation/SARPROZ/72Images/Brady_72img_TS_7.csv"
SarprozFile <- "d:/Geocomputation/SARPROZ/Brady_72img_TS_7.csv"
SarprozFile <- "d:/Geocomputation/SARPROZ/Brady_72img_Adjusted.csv"
PSInSAR_Jim <- read.csv2(SarprozFile, sep = ",", dec = ".", 
                     stringsAsFactors = F, strip.white = T)

PSInSAR <- subset(PSInSAR, select=c(LAT, LON, Vel01, COHER))
colnames(PSInSAR) <- c("Latitude", "Longitude", "Velocity", "Coherence")

PSInSAR_sf <- st_as_sf(PSInSAR, coords = c("Longitude", "Latitude"), 
                       agr = "constant", crs=4326)

Geo_ai_sf <- raster::extract(brady_ai_stack2, PSInSAR_sf, na.rm=TRUE)

# Loads and plots  thematic map of the state of Nevada
pacman::p_load(tigris)
nevada <- tigris::counties(state="NV", resolution="5m")
qtm(nevada[,"NAME"])

# n_cols <- length(names(PSInSAR))
# PSInSAR[1:3,(n_cols-2):n_cols]
# PSInSAR_coord <- PSInSAR[,c("X", "Y")]
# 
# date_cols <- grepl("X..[.]..[.]..", names(PSInSAR))
# PSInSAR_data <- PSInSAR[,date_cols]
# names(PSInSAR)[date_cols]
# 
# col_dates <- sub("X", "20", names(PSInSAR)[date_cols])
# col_dates <- as.numeric(format(as.Date(gsub("[.]", "-", col_dates)),"%Y%m%d"))
# names(PSInSAR_data)<- col_dates
# 
# names(PSInSAR)[1:10]

Velocities <- as.data.frame(PSInSAR[,c("X", "Y", "VEL")])
# Velocities$VEL <- as.numeric(Velocities$VEL)
rm(PSInSAR)

names(Velocities)
velocities_mean <- mean(Velocities$VEL)
velocities_sd <- sd(Velocities$VEL)
geo_factor <- 1
velocities_threshold <- velocities_mean - geo_factor*velocities_sd
nongeo_threshold <- velocities_mean + geo_factor*velocities_sd
Velocities$Geothermal <- ifelse(Velocities$VEL<velocities_threshold, 1, ifelse(Velocities$VEL>nongeo_threshold, 0, NA))

# s <- terra::rast("E:/HyMap_USA_2003_B_NV_HyVista/Brady_ref_L1C_mos.bil")
# s <- terra::rast(xmin=325209, xmax=336531, ymin=4395438, ymax=4412103, res=c(15,15), crs="epsg:32611")


s <- terra::rast(xmin=325209, xmax=336531, ymin=4395438, ymax=4412103, res=c(5,5), crs="epsg:32611")
s <- terra::rast(xmin=325209, xmax=(325209+1132*5), ymin=(4412103-1666*5), ymax=4412103, 
                 res=c(30,30), crs="epsg:32611")
# s <- terra::rast(deformation_layer)
v1 <- st_as_sf(Velocities, coords = c("X", "Y"), crs=crs(s))
v1 <- st_crop(v1, st_bbox(s))
r <- rasterize(vect(v1), terra::rast(s), field="VEL")
plot(r, col=heat.colors(20))
raster::spplot(raster::raster(r), main="Geothermal")



# rm(PSInSAR_coord, PSInSAR_data)
Velocities_pts <- st_multipoint(as.matrix(Velocities), dim="XYM")
head(Velocities)
fit_PSInSAR <- glm(VEL ~ X+Y, family = "gaussian", data = Velocities )
fit_PSInSAR

# s <- raster::stack("e:/HyMap_Brady/HyMap_Brady02.envi")
# # s <- s[[c(13,6,2)]]   # RGB
# s <- s[[c(6,2)]]
# plot(s[[1]])
# xm<-matrix(xFromCell(s,c(1:(raster::ncell(s)))),nrow=(raster::nrow(s)),byrow=TRUE)
# ym<-matrix(yFromCell(s,c(1:(raster::ncell(s)))),nrow=(raster::nrow(s)),byrow=TRUE)
# s[[1]] <- setValues(s[[1]], xm)
# s[[2]] <- setValues(s[[2]], ym)
# rm(xm, ym)
# names(s) <- c("X", "Y")
# pred <- raster::predict(s, model=fit_PSInSAR, type="response")
# plot(pred)
# rm(pred)

f_v <- Velocities
# library(raster)
sp::coordinates(f_v) <- ~ X + Y
sp::gridded(f_v) <- TRUE
sp::proj4string(f_v) <- crs(s)
sp::spplot(f_v)
f_v2 <- terra::rasterize(as.matrix(Velocities), s, field="VEL")
f_v2 <- f_v2[["VEL"]]

p_load(gstat)
sp <- as(s, "SpatialPoints")
k <- krige(VEL ~ x + y, f_v, sp)
Velocities_sf <- st_as_sf(Velocities, coords = c("X", "Y"), crs=32611) %>%
  cbind(st_coordinates(.))
summarize(Velocities_sf)
v2 <- st_crop(Velocities_sf, extent(s))
v_emp_OK <- gstat::variogram(
  VEL~1,
  as(v2, "Spatial") # switch from {sf} to {sp}
)
p_load(automap)
v_mod_OK <- automap::autofitVariogram(VEL~1, as(v2, "Spatial"))$var_model

grd_003_sf <- st_bbox(s) %>% st_as_sfc() %>% 
  st_make_grid(cellsize = c(3,3), what = "centers") %>%
  st_as_sf() %>%
  cbind(., st_coordinates(.))

grd_003_sp <- as(grd_003_sf, "Spatial")
gridded(grd_003_sp) <- TRUE
grd_003_sp <- as(grd_003_sp, "SpatialPixels")
