if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC)

# PSInSAR analysis
# Original file in: /store03/geodata/RawData/Mahmut/SAR/Results/PSI/Brady_and_Desert/2017-2019/72Images


doe_writeRaster <- function(x, filename, format="raster", overwrite=TRUE, bandorder="BSQ"){
  if(tools::file_ext(filename) != "grd") {
    filename <- tools::file_path_sans_ext(filename)
    filename <- paste(filename, ".grd", sep="")
  }
  f1<-raster::writeRaster(x=x, filename=filename, bandorder=bandorder, 
                          format=format, overwrite=overwrite)
  raster::hdr(f1, "ENVI")
  return(f1)
}


normalize_raster <- function(s){
  mnv <- raster::cellStats(s,'min', na.rm=TRUE)
  mxv <- raster::cellStats(s,'max', na.rm=TRUE)
  x <- (s - mnv) / (mxv - mnv)
  return(x)
}


som_stack_file = "d:/geoai_som/brady_som_output.gri"
som_stack_file = "d:/ThesisLayers/AI_Stacks/desert_ai_stack2.grd"
som_stack_file = "d:/ThesisLayers/AI_Stacks/brady_ai_stack2.grd"


roads <- raster::shapefile("E:/Jim+Erika/Built Up - projected.shp")
roads <- roads["id"]
som_stack <- raster::stack(som_stack_file)
deformation_layer <- raster::raster("D:/ThesisLayers/Deformation/Raster/Deformation.tif")

tcimf_minerals_layer <- raster::stack("D:/mm_final/HyMap/HyMapFull_tcimf")
tcimf_mineral_names <- names(minerals_layer)
tcimf_mineral_names <- gsub(".*\\.\\.(.*)\\.", "\\1", mineral_names)
cat("Minerals ok? ")
cat(paste0(
  all(tcimf_mineral_names == c("Chalcedony", "Kaolinite", "Gypsum", "Hematite", 
                             "Epsomite"))), 
  "\n")
tcimf_minerals_layer <- raster::projectRaster(tcimf_minerals_layer, som_stack)


osp_minerals_layer <- raster::stack("D:/mm_final/HyMap/HyMapFull_osp")
osp_mineral_names <- names(osp_minerals_layer)
osp_mineral_names <- gsub(".*\\.\\.(.*)\\.", "\\1", osp_mineral_names)
cat("Minerals ok? ")
cat(paste0(
  all(osp_mineral_names == c("Chalcedony", "Kaolinite", "Gypsum", "Hematite", 
                             "Epsomite"))), 
  "\n")
osp_minerals_layer <- raster::projectRaster(osp_minerals_layer, som_stack)
names(osp_minerals_layer) <- osp_mineral_names


sam_minerals_layer <- raster::stack("D:/mm_final/HyMap/HyMapFull_sam_rule")
sam_mineral_names <- names(minerals_layer)
sam_mineral_names <- gsub(".*\\.\\.(.*)\\.", "\\1", mineral_names)
cat("Minerals ok? ")
cat(paste0(
  all(sam_mineral_names == c("Chalcedony", "Kaolinite", "Gypsum", "Hematite", 
                         "Epsomite"))), 
  "\n")
sam_minerals_layer <- raster::projectRaster(sam_minerals_layer, som_stack)
sam_minerals_layer <- pi-sam_minerals_layer
minerals_layer[minerals_layer<0] <- 0
minerals_layer <- raster::scale(minerals_layer)
minerals_layer <- normalize_raster(minerals_layer)
names(minerals_layer) <- mineral_names

deformation_layer <- raster::projectRaster(deformation_layer, som_stack)

# for desert peak the threshold is -2 mm, for Brady it's ~ -0.87
geothermal_layer <- deformation_layer< (raster::cellStats(deformation_layer, mean)) # (-2) # (raster::cellStats(deformation_layer, mean))
names(geothermal_layer) <- "Geothermal"
d <- raster::scale(deformation_layer)
d <- normalize_raster(d)
names(d) <- "Deformation"

ai_stack <- raster::stack(geothermal_layer, 
                          som_stack[[c("Minerals", "Temperature", "Faults", "Slope")]],
                          d, #deformation_layer, 
                          minerals_layer)
m <- raster::mask(ai_stack, roads["id"], inverse = TRUE)
m <- normalize_raster(m)

doe_writeRaster(m, "d:/ThesisLayers/AI_Stacks/FullStacks/brady_full_stack")

m <- raster::stack("d:/ThesisLayers/AI_Stacks/FullStacks/brady_full_stack")

m2<-m[[c("Geothermal", "Temperature", "Faults", "Slope", "Chalcedony", 
         "Kaolinite", "Gypsum", "Hematite")]]
doe_writeRaster(m2, "d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack")
