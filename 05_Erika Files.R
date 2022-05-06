require("raster")
### START # Defining geographical extent from Areas of Interest ---------------
###

# Original Brady extent: (xmn=327499.1, xmx=329062.1, ymn = 4405906,
#                         ymx=4409320, res=c(3,3), crs="+init=epsg:32611")

extent_tall_brady <- raster::raster(
  xmn = 326440.1,
  xmx = 329062.1,
  ymn = 4405094,
  ymx = 4409320,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)

poly_extent_brady <- as(extent(extent_tall_brady), 'SpatialPolygons')
crs(poly_extent_brady) <- crs(extent_tall_brady)

extent_desert <- raster::raster(
  xmn = 330955.3,
  xmx = 335989.3,
  ymn = 4401681,
  ymx = 4406211,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)

poly_extent_desert <- as(extent(extent_desert), 'SpatialPolygons')
crs(poly_extent_desert) <- crs(extent_desert)

extent_hymap <- raster::raster(
  xmn = 325209,
  xmx = 336531,
  ymn = 4395438,
  ymx = 4412103,
  res = c(3, 3),
  crs = "+init=epsg:32611"
)
### END  # Defining geographical extent from Areas of Interest ----------------

selected_mineral <- "Chalcedony"
selected_method <- "mf"


a <- raster::stack(paste0("d:/erika/aster/Trial_12_", selected_method))
a <- raster::projectRaster(a, extent_tall_brady)
names(a) <- c("Chalcedony", "Hematite", "Epsomite", "Kaolinite", "Opal", "Gypsum")
s <- raster::stack(paste0("e:/mm_final/HyMap/HyMapFull_", selected_method))
s <- raster::projectRaster(s, extent_tall_brady)
names(s) <- c("Chalcedony", "Kaolinite", "Gypsum", "Hematite", "Epsomite") #, "OpalizedTuff", "KrattOpal")
# b <- raster::raster("e:/mm_final/HyMapBrady.gri")
a_ch <- a[[selected_mineral]]
s_ch <- s[[selected_mineral]]
# a_ch <- raster::projectRaster(a_ch, s_ch)
# mask_ch <- a_ch*0+1
# mask_ch[a_ch==0] <- NA
a_ch[a_ch==0] <- NA
s_ch[is.na(a_ch)] <- NA
spplot(stack(a_ch, s_ch))
th_s <- raster::quantile(s_ch, 0.99, na.rm=T)
th_a <- raster::quantile(a_ch, 0.99, na.rm=T)
s_ch[s_ch<th_s] <- NA
a_ch[a_ch<th_a] <- NA
th_a <- raster::quantile(a_ch, probs=0.99, na.rm=T)
th_s <- raster::quantile(s_ch, probs=0.99, na.rm=T)
a_ch[a_ch>th_a] <- th_a
s_ch[s_ch>th_s] <- th_s
raster::hist(a_ch)
raster::hist(s_ch)
spplot(a_ch, main=paste0("ASTER ", toupper(selected_method)))
spplot(s_ch, main=paste0("HyMap ", toupper(selected_method)))


