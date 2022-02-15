# ####
# remove.packages("spDataLarge")
# remove.packages("sp")
# remove.packages("sf")
# remove.packages("terra")
# ####


# install.packages("sf")
# install.packages("terra")
# install.packages("spData")
# install.packages("spDataLarge", repos = "https://nowosad.r-universe.dev")

if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, spData)
if(!require("spDataLarge")){
  install.packages("spDataLarge", repos = "https://nowosad.r-universe.dev")
  library("spDataLarge")
}







p_load(tmap)
data("lsl", package = "spDataLarge")
data("ta", package = "spDataLarge")
lsl_sf = st_as_sf(lsl, coords = c("x", "y"), crs = 32717)
hs = hillShade(ta$slope * pi / 180, terrain(ta$elev, opt = "aspect"))
rect = tmaptools::bb_poly(hs)
bbx = tmaptools::bb(hs, xlim = c(-0.02, 1), ylim = c(-0.02, 1), relative = TRUE)

tm_shape(hs, bbox = bbx) +
	tm_grid(col = "black", n.x = 1, n.y = 1, labels.inside.frame = FALSE,
	        labels.rot = c(0, 90)) +
	tm_raster(palette = gray(0:100 / 100), n = 100, legend.show = FALSE) +
	tm_shape(ta$elev) +
	tm_raster(alpha = 0.5, palette = terrain.colors(10), legend.show = FALSE) +
	tm_shape(lsl_sf) +
	tm_bubbles("lslpts", size = 0.2, palette = "-RdYlBu", title.col = "Landslide: ") +
  # tm_shape(sam) +
  # tm_bubbles(border.col = "gold", border.lwd = 2, alpha = 0, size = 0.5) +
  qtm(rect, fill = NULL) +
	tm_layout(outer.margins = c(0.04, 0.04, 0.02, 0.02), frame = FALSE) +
  tm_legend(bg.color = "white")
