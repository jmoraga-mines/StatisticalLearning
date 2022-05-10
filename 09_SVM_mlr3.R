############# START ### Once again, set up SP CV
if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC, mlr3spatiotempcv)
dir1 = "/store03/thesis/ai_stacks/MineralStacks"
dir2 = "d:/ThesisLayers/AI_Stacks/MineralStacks"
if(dir.exists(dir1)) {
  dir_base = dir1
} else if(dir.exists(dir2)){
  dir_base = dir2
} else {
  print("Error: data directories do not exist")
}
input_stack_file = file.path(dir_base, "brady_minerals_stack.gri")
m <- raster::stack(input_stack_file)
########       Spatioemporal sampling from here and on...
new_ai_data_xydf <- as.data.frame(m, xy=TRUE) # We will use coordinates for spatial sampling
new_ai_data_xydf <- na.omit(new_ai_data_xydf)
new_ai_data_xydf$Geothermal <- factor(ifelse(new_ai_data_xydf$Geothermal>=0.5, "Yes", "No"),
                                      levels = c("No", "Yes"))
# Creates a new taks to classify using spatio-temporal cross-validation
task_st_classif <- mlr3spatiotempcv::TaskClassifST$new(new_ai_data_xydf, 
                                                       target = "Geothermal",
                                                       positive = "Yes",
                                                       id = "geot_stcv",
                                                       extra_args = list(coords_as_features = FALSE,
                                                                         coordinate_names=c("x", "y"))
)
# Creates spatiotemporal resampling. 5 folds, repeated 2 times
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)
############# END  ### Once again, set up SP CV


# Support vector machine (SVM)
pacman::p_load(e1071)
svm_learner <- mlr3::lrn("classif.svm", 
                         cachesize = 1024, 
                         tolerance = 1, 
                         type = 'C-classification',
                         # gamma=0.2, # For kernels: polynomial, radial, sigmoid
                         cost=100,
                         kernel='linear') # SVM

print(svm_learner)
set.seed(42)
svm_sp <- mlr3::resample(task = task_st_classif,
                        learner = svm_learner,
                        #                        store_models = TRUE,
                        resampling = resampling_spcv)

svm_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
svm_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(rf_sp)

