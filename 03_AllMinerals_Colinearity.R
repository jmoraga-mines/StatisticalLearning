if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse)
pacman::p_load(fs)
pacman::p_load(DescTools)

mineral_thresholds <- c(0.12, 0.045, 0.04, 0.055, 0.02, 0.06, 0.025)

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
circular_weight <- function(radius, crs=32611, fillNA=FALSE, cell_value=NA) {
  empty_raster <- raster::raster(ncols=21, nrows=21, xmn=0, res=c(1, 1), crs=crs)
  Weight_matrix <- raster::focalWeight(empty_raster, radius, "circle", fillNA = fillNA)
  if(!is.na(cell_value)){
    center_value <- as.matrix(Weight_matrix)[radius+1, radius+1]
    Weight_matrix <- Weight_matrix*cell_value/center_value
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


# DON'T RUN
# raster::crs(original_envi_stack) <- 32611
# raster::crs(original_envi_stack) <- raster::crs("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs")
# original_envi_stack <- original_envi_stack[[1:7]]
# minerals_stack <- rast(minerals_file)
minerals_stack <- read_envi_file(original_file_name, minerals_file_name)
print(minerals_stack)


# DON'T RUN
# terra::crs(minerals_stack) <- "epsg:32611"
# terra::res(minerals_stack) <- c(3,3)
# terra::ext(minerals_stack) <- terra::ext(original_envi_stack)
# 
# names(minerals_stack) <- sub("Target [(]TCIMF Score [(]([A-Za-z_ ]+)[)][)]", # Regular Expression
#                              "\\1",  # Replace by the first match within parentheses 
#                              names(minerals_stack))
# input_stack_file = "d:/geoai_som/brady_som_output.gri"
# som_output_stack <- terra::rast(input_stack_file)

# minerals_projected <- terra::project(minerals_stack, som_output_stack, method="near")
projected_minerals_stack <- raster::projectRaster(minerals_stack, som_output_stack)
projected_minerals_stack <- raster::stack(projected_minerals_stack)

# values(original_envi_stack) <- values(minerals_stack)
# names(original_envi_stack) <- sub("[A-Za-z.]+[.]{2}([A-Za-z_ ]+)[.]{1,2}$", # Regular Expression
#                             "\\1",  # Replace by the first match within parentheses 
#                             names(minerals_stack))

minerals_stack_df <- as.data.frame(projected_minerals_stack)

# DO NOT RUN
# Unselects unneeded minerals
# minerals_stack_df <- <- minerals_stack_df %>% select(!c(Opalized_Tuff, KrattOpal))
print(head(minerals_stack_df))

cov(minerals_stack_df)
opals <- projected_minerals_stack[[c("Chalcedony", "Opalized_Tuff","KrattOpal")]]
opals <- setMinMax(opals)
opals_threshold <- opals-c(0.12, 0.06, 0.025)
opals_threshold[opals_threshold<0] <- 0
opals_scaled <- normalize_raster(opals_threshold)
opals_df <- as.data.frame(opals_scaled)
opals_winsor <- as.data.frame(lapply(opals_df, winsor_right))
plot(opals_scaled)
opals_winsor_normal <- as.data.frame(lapply(opals_winsor, normalize)) 
cor(opals_df)
cor(opals_winsor)


W1 <- circular_weight(30, fillNA = T)
s1_focal <- terra::focal(terra::rast(projected_minerals_stack), W1)
#names(s1_focal) <- names(minerals_stack)
s2_focal <- raster::stack(s1_focal)
s3_focal <- raster::stack(som_output_stack, s2_focal)
s3_focal <- s3_focal[[c("Geothermal", "Minerals", "Temperature", "Faults",
                        "Chalcedony", "Kaolinite", "Gypsum", "Hematite",
                        "Epsomite")]]
s3_focal_df <- as.data.frame((s3_focal))
rm(r, W1, s1_focal, s2_focal)

# Discard tuples with NA values
s3_focal_cor <- cor(s3_focal_df, use = 'complete')
s3_focal_cov <- cov(s3_focal_df, use = 'complete')

################## Ok ,let's add the layers in terra
som_output_stack_terra <- project(rast(minerals_stack), som_output_stack)
som_output_stack_terra <- c(som_output_stack, som_output_stack_terra)

model_basic <- glm(Geothermal ~ .-Opalized_Tuff-KrattOpal-Subsidence-Uplift, 
                   family = binomial, data = s3_df)
s3_cor <- cor(s3_df)
s3_cov <- cov(s3_df)
write.csv(as.data.frame(s3_cor), "all_layer_correlation.csv")
write.csv(as.data.frame(s3_cov), "all_layer_covariance.csv")
vif_original <- car::vif(model_basic)

s3_df <- as.data.frame(som_output_stack_terra)
s3_focal_df <- as.data.frame(s3_focal)

model_focal_basic <- glm(Geothermal ~ .-Opalized_Tuff-KrattOpal-Subsidence-Uplift,
                         family = binomial, data = s3_focal_df)
s3_focal_cor <- cor(s3_focal_df)
s3_focal_cov <- cov(s3_focal_df)
write.csv(as.data.frame(s3_cor), "all_layer_focal_correlation.csv")
write.csv(as.data.frame(s3_cov), "all_layer_focal_covariance.csv")
vif_focal <- car::vif(model_focal_basic)
model_TFKGE <- glm(Geothermal ~ Temperature+Faults+Kaolinite+Gypsum+Epsomite,
                   family = binomial, data = s3_df)
model_TFCHE <- glm(Geothermal ~ Temperature+Faults+Chalcedony+Hematite+Epsomite,
                   family = binomial, data = s3_df)
model_TFCG  <- glm(Geothermal ~ Temperature+Faults+Chalcedony+Gypsum,
                   family = binomial, data = s3_df)
vif_TFKGE <- car::vif(model_TFKGE)
vif_TFCHE <- car::vif(model_TFCHE)
vif_TFCG  <- car::vif(model_TFCG)

