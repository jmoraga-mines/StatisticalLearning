if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, dplyr, caret)

input_stack_file = "d:/geoai_som/brady_som_output.gri"
s <- terra::rast(input_stack_file)
cat("Stack loaded:")
print(s)
plot(s[[1]])
stack_sf <- as.data.frame(s, xy = TRUE)
stack_sf$Geothermal <- sapply(stack_sf$Geothermal, 
                              function(x) (if (x==0) return("FALSE") 
                                           else return("TRUE")))
stack_sf$Geothermal <- factor(stack_sf$Geothermal, levels = c("FALSE", "TRUE"))
pacman::p_load(caret, tidyverse)
set.seed(123)
training.samples <- stack_sf$Geothermal %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- stack_sf[training.samples, ]
test.data <- stack_sf[-training.samples, ]
train.data <- stack_sf
# Build the models
### Build glm model with 3 inputs
model1 <- glm(Geothermal ~ Minerals + Temperature + Faults, 
              family = binomial, data = train.data)
cat("AUROC: model 1\n")
print(pROC::auc(pROC::roc(train.data$Geothermal, fitted(model1))))
### Build glm model with 3 inputs and coordinates
model2 <- glm(Geothermal ~ x + y + Minerals + Temperature + Faults, 
              family = binomial, data = train.data)
cat("AUROC: model 2 (includes x, y coordinates)\n")
print(pROC::auc(pROC::roc(train.data$Geothermal, fitted(model2))))

# Make predictions
predictions1 <- model1 %>% predict(test.data)
predictions2 <- model2 %>% predict(test.data)
p1 <- sapply(predictions1,
             function(x) (if (x<0.5) return("FALSE")
                          else return("TRUE")))
p1 <- factor(p1, levels =  c("FALSE", "TRUE"))
p2 <- sapply(predictions2,
             function(x) (if (x<0.5) return("FALSE")
                          else return("TRUE")))
p2 <- factor(p2, levels = c("FALSE", "TRUE"))
# Model performance
pacman::p_load(car)
cat("Model 1\n")
print(model1)
data.frame(
  RMSE = RMSE(predictions1, test.data$Geothermal),
  R2 = R2(predictions1, test.data$Geothermal)
)
cat("Model 1 VIF\n")
car::vif(model1)

cat("\n\nModel 2\n")
print(model2)
print(data.frame(
  RMSE = RMSE(predictions2, test.data$Geothermal),
  R2 = R2(predictions2, test.data$Geothermal)
))
cat("Model 2 VIF\n")
print(car::vif(model2))

############### create raster from model 1
t <- s[[c("Minerals", "Temperature", "Faults")]]
pred1 <- raster::predict(t, model = model1, type="response")
plot(c(pred1,s[[1]]), main=c("Model 1 prediction", "Ground truth"))
plot(c(pred1>0.5,s[[1]]), main=c("Model 1 pred [0,1]", "Ground truth"))

############### create raster from model 2
t <- s[[c("Minerals", "Temperature", "Faults", "Subsidence", "Uplift")]]
xm<-matrix(xFromCell(t,c(1:(raster::ncell(t)))),nrow=(raster::nrow(t)),byrow=TRUE)
ym<-matrix(yFromCell(t,c(1:(raster::ncell(t)))),nrow=(raster::nrow(t)),byrow=TRUE)
t[[4]] <- setValues(t[[4]], xm)
t[[5]] <- setValues(t[[5]], ym)
names(t) <- c("Minerals", "Temperature", "Faults", "x", "y")
pred2 <- raster::predict(t, model = model2, type="response")
plot(c(pred2,s[[1]]), main=c("Model 2 prediction", "Ground truth"))
plot(c(pred2>0.5,s[[1]]), main=c("Model 2 pred [0,1]", "Ground truth"))

plot(c(pred1, pred2), main=c("Model 1 prediction", "Model 2 prediction"))
# plot(c(pred1, pred2,s[[1]]), main=c("Model 1 pred", "Model 2 pred", "Ground truth"))


################################################################################
# stack_sf <- st_as_sf(stack_sf, coords = c("x", "y"), crs=crs(s))
# stack_sf <- stack_sf %>% 
#   select(c("Geothermal", "Minerals", "Temperature", "Faults"))
# head(stack_sf)
# summary(stack_sf)

################################################################################

pacman::p_load(mlr)
stack_coords <- stack_sf[,c("x", "y")]
stack_data <- dplyr::select(stack_sf, -x, -y)
# stack_data$Geothermal <- factor(stack_data$Geothermal, levels = 0:1)
# create task
task <- mlr::makeClassifTask(data = stack_data, target = "Geothermal",
                             positive = "TRUE", coordinates = stack_coords)
binomial_learner <- makeLearner(cl = "classif.binomial",
                                link = "logit",
                                predict.type = "prob",
                                fix.factors.prediction = TRUE)
spatial_resample <- makeResampleDesc(method = "SpRepCV", # Spatial repeated cross-validation
                                     folds = 5, 
                                     reps = 10) # Was 100 reps
set.seed(012345)
# set.seed(as.numeric(Sys.time()))
sp_cv = mlr::resample(learner = binomial_learner, task = task,
                      resampling = spatial_resample, 
                      measures = list(mlr::auc, # Measures AUC 
                                      mlr::acc, # Measures accuracy
                                      mlr::bac, # Measures binary accuracy
                                      mlr::tp, 
                                      mlr::tn,
                                      mlr::fp,
                                      mlr::fn) 
                      )
save(sp_cv, file = "sp_cv_50_iter.RData")
load(file = "sp_cv_50_iter.RData")

summary(sp_cv$measures.test$auc)
mean(sp_cv$measures.test$auc)
summary(sp_cv$measures.test$acc)
mean(sp_cv$measures.test$acc)
