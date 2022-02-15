if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse)
pacman::p_load(fs)
pacman::p_load(DescTools)

mineral_thresholds <- c(0.12, 0.045, 0.04, 0.055, 0.02, 0.06, 0.025)
layers_to_keep <- c("Geothermal", "Temperature", "Faults", "Chalcedony", 
                    "Kaolinite", "Gypsum", "Hematite", "Epsomite")

### Spectral analysis files from ENVI have missing CRS and extension
### information, this function fixes this error by using the original
### file analyzed to get CRS, resolution and extension data
read_envi_file <- function(original_file_name, analysis_file_name){
  if(!fs::file_exists(original_file_name)){
    cat(paste0("ERROR: Original file '", original_file_name, "' does not exist"))
    return(NULL)
  }
  if(!fs::file_exists(analysis_file_name)){
    cat(paste0("ERROR: Analysis file '", analysis_file_name, "' does not exist"))
    return(NULL)
  }
  analysis_stack <- raster::stack(analysis_file_name)
  new_names <- sub("[A-Za-z.]+[.]{2}([A-Za-z_ ]+)[.]{1,2}$", # Regular Expression
                                           "\\1",  # Replace by the first match within parentheses 
                                           names(analysis_stack))

  # In case the band names start with "Band.<num>" where <num> is the Band number
  if(grepl("Band[.][0-9]+", new_names[1])==TRUE)
    new_names <- sub("(Band[.][0-9]+)[.].*", "\\1", new_names)

  names(analysis_stack) <- new_names
  if(is.na(crs(analysis_stack))){
    cat("Analysis file has missing information, using original file to fix\n")
    original_stack <- raster::raster(original_file_name)
    if(is.na(crs(original_stack))){
      cat("ERROR: Original file has missing CRS information\n")
      return(NULL)
    }
    if(!((nrow(analysis_stack)==nrow(original_stack)) &&
         (ncol(analysis_stack)==ncol(original_stack)))){
      cat("ERROR: Original file and analysis file have different dimensions\n")
      return(NULL)
    }
    s <- raster::stack()
    for (layer in 1:(raster::nlayers(analysis_stack))){
      raster::values(original_stack) <- raster::values(analysis_stack[[layer]])
      names(original_stack) <- names(analysis_stack[[layer]])
      s <- raster::stack(s, original_stack)
    }
    analysis_stack <- s
    rm(s)
  }
  return(analysis_stack)
}

##########
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

winsor_raster <- function(s, probs=0.95){
  winsor_max <- raster::calc(s,'quantile', probs = probs, na.rm=TRUE)
  s <- raster::calc(s, fun = function(x, na.rm) pmin(x,winsor_max))
  return(s)
}

normalize_raster <- function(s){
  mnv <- raster::cellStats(s,'min', na.rm=TRUE)
  mxv <- raster::cellStats(s,'max', na.rm=TRUE)
  x <- (s - mnv) / (mxv - mnv)
  return(x)
}

normalize <- function(s){
  min_val <- min(s, na.rm = TRUE)
  max_val <- max(s, na.rm = TRUE)
  x <- (s - min_val) / (max_val - min_val)
  return(x)
}

# for vectors/lists
winsor_right <- function(x, probs=0.95, rm.lt_zero=TRUE, rm.na=TRUE){
  y = x
  if(rm.lt_zero)
    y <- y[y>0]
  if(rm.na)
    y <- y[!is.na(y)]
  winsor_max <- quantile(y, probs = probs)
  x[x>winsor_max] <- winsor_max
  return(x)
}

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
# Create square weight matrix
square_weight <- function(width, fillNA=FALSE, cell_value=NA, 
                            zero_center=FALSE) {
  if((width%%2)==0){
    cat("weights width must be an odd number\n")
    return(NULL)
  }
  Weight_matrix <- matrix(1,nrow=width, ncol=width)
  if(zero_center){
    Weight_matrix[floor(width/2)+1, floor(width/2)+1] = 0
  }
  return(Weight_matrix)
}

###### Input file, has SOM results for Geothermal (AI input)
som_stack_file = "d:/geoai_som/brady_som_output.gri"
som_output_stack <- raster::stack(som_stack_file)
som_output_stack <- som_output_stack[[c("Geothermal", "Minerals", "Temperature", "Faults")]]
cat("Cropping SOM file to only keep relevant layers (eliminate SUbsidence and Uplift)\n")
print(names(som_output_stack))

###### Get mineral analysis results from TCIMF
minerals_directory <- 'd:/mm_final/Brady'
minerals_file_name <- path(minerals_directory, 'HyMapBrady_tcimf')

original_file_name <- "d:/mm_final/HyMapBrady.gri"
original_envi_stack <- raster::stack(original_file_name)

# Read ENVI output for target spectra analysis
minerals_stack <- read_envi_file(original_file_name, minerals_file_name)
print(minerals_stack)

# Project (and crop), to area of interest
projected_minerals_stack <- raster::projectRaster(minerals_stack, som_output_stack)
projected_minerals_stack <- raster::stack(projected_minerals_stack)
new_ai_data <- raster::stack(som_output_stack, projected_minerals_stack)
new_ai_data <- new_ai_data[[layers_to_keep]]
minerals_stack_df <- as.data.frame(projected_minerals_stack)
new_ai_data_df <- as.data.frame(new_ai_data)

print(head(minerals_stack_df))
print(head(new_ai_data_df))
plot(new_ai_data)

# Discard tuples with NA values
minerals_stack_df_cor <- cor(minerals_stack_df, use = 'complete')
minerals_stack_df_cov <- cov(minerals_stack_df, use = 'complete')

################## Ok ,let's add the layers in terra
model_TF <- glm(Geothermal ~ Temperature + Faults, family = binomial, data = new_ai_data_df)
new_ai_data_df_cor <- cor(new_ai_data_df, use = 'complete')
new_ai_data_df_cov <- cov(new_ai_data_df, use = 'complete')

write.csv(as.data.frame(new_ai_data_df_cor), "all_layer_correlation.csv")
write.csv(as.data.frame(new_ai_data_df_cov), "all_layer_covariance.csv")
vif_original <- car::vif(model_TF)
cat("VIF for Generalized Linear Model:\n")
print(vif_original)

s3_df <- as.data.frame(new_ai_data)
model_TFKGE <- glm(Geothermal ~ Temperature+Faults+Kaolinite+Gypsum+Epsomite,
                   family = binomial, data = s3_df)
model_TFCHE <- glm(Geothermal ~ Temperature+Faults+Chalcedony+Hematite+Epsomite,
                   family = binomial, data = s3_df)
model_TFCG  <- glm(Geothermal ~ Temperature+Faults+Chalcedony+Gypsum,
                   family = binomial, data = s3_df)

vif_TFKGE <- car::vif(model_TFKGE)
cat("VIF for TFKGE Model:\n")
print(vif_TFKGE)
vif_TFCHE <- car::vif(model_TFCHE)
cat("VIF for TFCHE Model:\n")
print(vif_TFCHE)
vif_TFCG  <- car::vif(model_TFCG)
cat("VIF for TFCG Model:\n")
print(vif_TFCG)

###################
pacman::p_load(caret)

ROC_basic <- pROC::roc(model_TF$y, fitted(model_TF))
ROC_TFCG  <- pROC::roc(model_TFCG$y, fitted(model_TFCG))
ROC_TFKGE <- pROC::roc(model_TFKGE$y, fitted(model_TFKGE))
ROC_TFCHE <- pROC::roc(model_TFCHE$y, fitted(model_TFCHE))
AUROC_basic <- pROC::auc(ROC_basic)*100
AUROC_TFCG  <- pROC::auc(ROC_TFCG)*100
AUROC_TFKGE <- pROC::auc(ROC_TFKGE)*100
AUROC_TFCHE <- pROC::auc(ROC_TFCHE)*100

TF_glm_raster    <- raster::predict(new_ai_data, model_TF)
TFCG_glm_raster  <- raster::predict(new_ai_data, model_TFCG)
TFKGE_glm_raster <- raster::predict(new_ai_data, model_TFKGE)
TFCHE_glm_raster <- raster::predict(new_ai_data, model_TFCHE)
names(TF_glm_raster)    <- "TF_glm_raster"
names(TFCG_glm_raster)  <- "TFCG_glm_raster"
names(TFKGE_glm_raster) <- "TFKGE_glm_raster"
names(TFCHE_glm_raster) <- "TFCHE_glm_raster"
glm_rasters <- raster::stack(TF_glm_raster, TFCG_glm_raster, TFKGE_glm_raster, TFCHE_glm_raster)
plot(glm_rasters>0.5)
plot(new_ai_data[["Geothermal"]])

## Plotting results
pacman::p_load(tmap, tmaptools)
tmap_mode("plot")
### TIGER/Line and cartographic boundary shapefiles from the US Census Bureau 
pacman::p_load(tigris)

### Download Nevada Shapefile
nevada <- tigris::counties(state="NV", resolution="100m")
### Download Churchill County, NV shapefile
Churchill <- tigris::county_subdivisions("NV", county = "Churchill", resolution="10m")
# DON'T RUN
# Plots Nevada and counties by name
# tm_shape(nevada[,"NAME"])+tm_polygons(col="NAME")+tm_legend(legend.outside=T)

### Select county subdivision from Churchill, NV
CarsonSink <- Churchill %>% filter(NAME=="Carson Sink") # %>% select(NAME) 

# tm_shape(new_ai_data[["Geothermal"]])+tm_raster(style="cont", palette = "Greens")+tm_legend(show=T, legend.outside=T)
# Plots Nevada, highlighting Churchill County and subdivisions
nevada_counties <- tm_shape(nevada[,"NAME"])+tm_polygons(col="white", legend.show = F)+
  tm_shape(Churchill[,"NAME"])+
  tm_polygons(lwd=1, col="NAME", title = "Churchill County", border.col = "black")+
  tm_legend(show=T, legend.outside=T)+
  tm_shape(st_as_sfc(st_bbox(new_ai_data)))+
  tm_polygons(lwd=2, col="red", alpha = 0)
# Plots Carson Sink and the location of the Brady raster
tm_shape(CarsonSink$geometry)+tm_polygons(col = "grey")+
  tm_shape(new_ai_data[["Geothermal"]])+
  tm_raster(style="cat", palette = "Greens", labels=c("False", "True"))+
  tm_legend(show=T, legend.outside=T, title = 'Carson Sink') +
  tm_grid()+
  tm_compass(type="arrow")+
  tm_scale_bar()

# tm_shape(new_ai_data[["Geothermal"]])+
#   tm_raster(style="cat", palette = c("#DDDDDD", "#44AA44"), labels = c("False", "True"))+
#   tm_legend(show=T, legend.outside=T, title = "Brady")+
#   tm_compass(type="arrow", text.size = 0.8, text.color = "Black", 
#              position = c(0.005, .90)) +
#   tm_scale_bar(position = c(0.6, 0.005), text.size=1.2)

tm_shape(new_ai_data[["Geothermal"]])+
  tm_raster(style="cat", palette = c("#DDDDDD", "#44AA44"), labels = c("False", "True"))+
  tm_grid() +
  tm_legend(show=T, legend.outside=T, title = "Brady")+
  tm_compass(type="arrow", text.size = 0.8, text.color = "Black", 
             position = c("left", "top")) +
  tm_scale_bar(position = c("right", "bottom"), text.size=1.2) +
  tm_layout(inner.margins = 0)



pacman::p_load(OpenStreetMap)
osm_map <- tmaptools::read_osm(st_bbox(Churchill))
new_ai_data2 <- new_ai_data
new_ai_data2[new_ai_data<=0] <- NA

tmap_mode("plot")
tm_shape(osm_map) + tm_rgb()+ 
  tm_shape(new_ai_data2[["Chalcedony"]])+
  tm_raster(alpha=.4, style="quantile", palette = "YlOrRd")+
  tm_legend(show=T, legend.outside=T)+
  tm_compass(type="arrow", text.size = 0.8, text.color = "Black", position = c(0.005, .85)) +
  tm_scale_bar(position = c(0.6,0.05, size=1.2))

# tm_raster(new_ai_data)
moran_chalcedony <- list()
for (my_radius in 11:40){
  w1 <- circular_weight(my_radius,cell_value = 1, zero_center = TRUE)
  moran_chalcedony[my_radius] <- raster::Moran(new_ai_data[["Chalcedony"]], w1)
  print(moran_chalcedony[my_radius])
}
