if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, tidyverse, mlr3verse, nnet)

dir1 = "/store03/geodata/MachineLearning/doe-som"
dir2 = "d:/geoai_som"
if(dir.exists(dir1)) {
  dir_base = dir1
} else if(dir.exists(dir2)){
  dir_base = dir2
} else {
  print("Error: data directories do not exist")
}

input_stack_file = file.path(dir_base, "brady_som_output.gri")
test_stack_file = file.path(dir_base, "desert_som_output.gri")

s <- raster::stack(input_stack_file)
cat("Stack loaded:")
print(s)
s <- s[[c("Geothermal", "Minerals", "Temperature", "Faults")]]

stack_sf <- as.data.frame(s, xy = TRUE)
stack_sf <- na.omit(stack_sf)
stack_sf$Geothermal <- factor(if_else(stack_sf$Geothermal>=0.5, "true", "false"),
                              levels = c("true", "false"))
stack_sf <- st_as_sf(stack_sf, coords = c("x", "y"), crs=32611)
set.seed(123)
task_classification <- mlr3spatiotempcv::TaskClassifST$new(stack_sf, target = "Geothermal",
                                                           positive = "true", 
                                                           id = "geothermal_task",
                                                           extra_args = list(coords_as_features = FALSE, crs=crs(s))
)
resampling_sp <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)
nnet_learner <- mlr3::lrn("classif.nnet", size = 1) # Neural network, with hidden neurons

print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_classification,
                        learner = nnet_learner,
                        store_models = TRUE,
                        resampling = resampling_sp)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))


test_stack <- raster::stack(test_stack_file)
cat("Stack loaded:")
print(test_stack)
test_stack <- test_stack[[c("Geothermal", "Minerals", "Temperature", "Faults")]]

test_df <- as.data.frame(test_stack, xy = FALSE)
test_df <- na.omit(test_df)
test_df$Geothermal <- factor(if_else(test_df$Geothermal>=0.5, "true", "false"),
                          levels = c("true", "false"))

for (i in seq(1,10)){
  cat(paste0("Model nnet: ", i, "\n"))
  test_prediction <- predict(nn_sp$learners[[i]]$model, newdata = test_df, predict_type="prob") 
  test_prediction <- factor(if_else(test_prediction<0.5, "true", "false"),
                            levels = c("true", "false"))
  print(caret::confusionMatrix(data= test_prediction, reference= test_df$Geothermal))
}

