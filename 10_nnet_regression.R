if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC)
brady_ai_stack <- raster::stack("d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack")
brady_full_stack <- raster::stack("d:/ThesisLayers/AI_Stacks/FullStacks/brady_full_stack")
brady_ai_stack$Geothermal <- brady_full_stack$Deformation
rm(brady_full_stack)

new_ai_data_xydf <- as.data.frame(brady_ai_stack, xy=TRUE) # We will use coordinates for spatial sampling
new_ai_data_xydf <- na.omit(new_ai_data_xydf)
task_st_regr <- mlr3spatiotempcv::TaskRegrST$new(new_ai_data_xydf,
                                                 target = "Geothermal",
                                                 id = "geo_reg_stcv",
                                                 extra_args = list(coords_as_features = FALSE,
                                                                   coordinate_names=c("x", "y"))
                                                 )
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)

nnet_learner <- mlr3::lrn("classif.nnet", size = 7) # Neural network, with hidden neurons

print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_st_classif,
                        learner = nnet_learner,
                        #                        store_models = TRUE,
                        resampling = resampling_spcv)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

                                                       
