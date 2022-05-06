if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, tidyverse, mlr3verse)

input_som_file = "d:/geoai_som/brady_som_output.gri" # "d:/ThesisLayers/AI_Stacks/brady_ai_stack2.grd"
input_stack_file = "d:/ThesisLayers/AI_Stacks/brady_ai_stack2.grd"
input_stack_file = "d:/geoai_som/brady_som_output.gri"
# b <- raster::raster(input_som_file)
s <- raster::stack(input_stack_file)
# s <- raster::projectRaster(s, b)
cat("Stack loaded:")
print(s)
s <- s[[c("Geothermal", "Minerals", "Temperature", "Faults")]]
#i <- is.na(s)
#s[i] <- 0
stack_sf <- as.data.frame(s, xy = TRUE)
stack_sf <- na.omit(stack_sf)
stack_sf$Geothermal <- factor(if_else(stack_sf$Geothermal==1, "true", "false"),
                              levels = c("true", "false"))
#stack_sf$row_id <- 1:nrow(stack_sf)
#geothermal_backend <- as_data_backend(stack_sf, primary_key = "row_id")
stack_sf <- st_as_sf(stack_sf, coords = c("x", "y"), crs=32611)
set.seed(123)
task_classification <- mlr3spatiotempcv::TaskClassifST$new(stack_sf, target = "Geothermal",
                                                           positive = "true", 
                                                           id = "geothermal_task",
                                                           extra_args = list(coords_as_features = FALSE)
                                                           )
binary_learner <- mlr3::lrn("classif.log_reg", id = "binary") # Logistic regression
resampling_sp <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)
set.seed(42)
rr_sp <- mlr3::resample(task = task_classification,
                        learner = binary_learner,
                        store_models = TRUE,
                        resampling = resampling_sp)
print(rr_sp)
# rr_sp$score()
rr_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
rr_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rr_sp$aggregate(msr("classif.ce"))
rr_sp$aggregate(msr("classif.acc"))
rr_sp$aggregate(msr("classif.bacc"))
rr_sp$aggregate(msr("classif.fbeta"))
rr_sp$aggregate(msr("classif.sensitivity"))
rr_sp$aggregate(msr("classif.specificity"))

rr_sp$aggregate(msr("classif.ce", id = "ce.train", predict_sets = "train"))
rr_sp$aggregate(msr("classif.ce", id = "ce.test", predict_sets = "test"))
for(i in 1:10){
  cat(paste0("Model:", i, "\n"))
  print(rr_sp$learners[[i]]$model$score(msr("classif.bacc")))
  print(rr_sp$learners[[i]]$model$score(msr("classif.fbeta")))
  print(rr_sp$learners[[i]]$model$coefficients)
  cat(paste0("AIC:", rr_sp$learners[[i]]$model$aic, "\n"))
}

pred=rr_sp$prediction()
print(pred$confusion)
print("Accuracy:")
print(pred$score(msr("classif.acc")))
print("Sensitivity:")
print(pred$score(msr("classif.sensitivity")))
print("Specificity:")
print(pred$score(msr("classif.specificity")))
print("Balanced Accuracy:")
print(pred$score(msr("classif.bacc")))
print("f-1:")
print(pred$score(msr("classif.fbeta")))
autoplot(pred, type = "roc")


partitions_plot <- autoplot(resampling_sp, task_classification, 
                            fold_id = c(1:5), size = 0.3, 
                            crs=st_crs(32611))
plot(partitions_plot)


############## Don't run below


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
stack_sf <- as.data.frame(s, xy = FALSE)
stack_sf <- na.omit(stack_sf)
stack_sf$Geothermal <- factor(if_else(stack_sf$Geothermal==1, "true", "false"),
                              levels = c("true", "false"))
resampling_sp2 <- mlr3::rsmp('repeated_cv', folds = 5, repeats = 2)
task_classification2 <- mlr3::as_task_classif(stack_sf,
                                            target = "Geothermal",
                                            positive = "true",
                                            id = "geothermal_task")
resampling_sp2 <- mlr3::rsmp('repeated_cv', folds = 5, repeats = 2)
logistic_learner <- mlr3::lrn("classif.log_reg", id = "binary") # Logistic regression
print(logistic_learner)


rr_sp2 <- mlr3::resample(task = task_classification2,
                        learner = logistic_learner,
                        store_models = TRUE,
                        resampling = resampling_sp2)
rr_sp2$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
rr_sp2$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))


train_set <- sample(task_classification2$nrow, 0.8 * task_classification2$nrow)
test_set  <- setdiff(seq_len(task_classification2$nrow), train_set)
logistic_learner$train(task_classification2, row_ids = train_set) # NULL) #, row_ids = train_set)
prediction <- logistic_learner$predict(task_classification2, row_ids = test_set) # row_ids = NULL) #row_ids = test_set)
print(prediction$confusion)
print(prediction$score(msr("classif.acc")))
print(prediction$score(msr("classif.bacc")))
print(prediction$score(msr("classif.fbeta")))
print(prediction$score(msr("classif.sensitivity")))
print(prediction$score(msr("classif.specificity")))
print(prediction$score(msr("classif.auc")))
autoplot(prediction, type = "roc")
