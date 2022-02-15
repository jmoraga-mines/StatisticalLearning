if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, tidyverse, mlr3verse)

input_stack_file = "d:/geoai_som/brady_som_output.gri"
s <- terra::rast(input_stack_file)
cat("Stack loaded:")
print(s)
s <- s[[c("Geothermal", "Minerals", "Temperature", "Faults")]]
stack_sf <- as.data.frame(s, xy = TRUE)
stack_sf$Geothermal <- factor(if_else(stack_sf$Geothermal==1, "true", "false"),
                              levels = c("false", "true"))
# stack_sf$row_id <- 1:nrow(stack_sf)
geothermal_backend <- as_data_backend(stack_sf, primary_key = "row_id")
stack_sf <- st_as_sf(stack_sf, coords = c("x", "y"), crs=32611)
task_classification <- mlr3spatiotempcv::TaskClassifST$new(stack_sf, target = "Geothermal",
                                                           positive = "true", 
                                                           id = "geothermal_task",
                                                           extra_args = list(coords_as_features = FALSE)
                                                           )
binary_learner <- mlr3::lrn("classif.log_reg", id = "binary") # Logistic regression
resampling_sp <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)
rr_sp <- mlr3::resample(task = task_classification,
                        learner = binary_learner,
                        resampling = resampling_sp)

partitions_plot <- autoplot(resampling_sp, task_classification, 
                            fold_id = c(1:5), size = 0.3, 
                            crs=st_crs(32611))
plot(partitions_plot)


### Errors below
# autoplot(resampling_sp, task, fold_id = c(1:4), size = 0.7) *
#   ggplot2::scale_y_continuous(breaks = seq(-3.97, -4, -0.01)) *
#   ggplot2::scale_x_continuous(breaks = seq(-79.06, -79.08, -0.01))

# task_classification <- mlr3spatiotempcv::TaskClassifST$new(stack_sf, target = "Geothermal",
#                                          positive = "true", 
#                                          id = "geothermal_task",
#                                          extra_args = list(coords_as_features = FALSE, 
#                                                            crs = crs(s), 
#                                                            coordinate_names = c("x", "y")))
# task_classification <- mlr3::as_task_classif(stack_sf,
#                                              target = "Geothermal",
#                                              positive = "true",
#                                              id = "geothermal_task")
binary_learner <- mlr3::lrn("classif.log_reg", id = "binary") # Logistic regression
print(binary_learner)

train_set <- sample(task_classification$nrow, 0.8 * task_classification$nrow)
test_set  <- setdiff(seq_len(task_classification$nrow), train_set)
binary_learner$train(task_classification, row_ids = train_set)
prediction <- binary_learner$predict(task_classification, row_ids = test_set)
print(prediction$confusion)
print(prediction$score(msr("classif.acc")))
print(prediction$score(msr("classif.bacc")))
autoplot(prediction, type = "roc")
