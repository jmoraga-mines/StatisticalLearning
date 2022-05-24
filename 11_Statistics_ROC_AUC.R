library(sf)
b <- raster::raster("d:/ThesisLayers/Normalized/HyMap_brady_ai_stack")
d <- raster::raster("d:/ThesisLayers/Normalized/HyMap_desert_ai_stack")
p <- raster::raster("d:/ThesisLayers/Normalized/prediction/HyMapD_27x7.13.grd")
df_results <- cbind(values(p), values(d))
colnames(df_results) <- c("y", "y_pred")
df_results <- na.omit(df_results)
df_results <- as.data.frame(df_results)
df_results$y_pred = factor(df_results$y_pred,levels= c("0", "1"), labels=c("No", "Yes"))
df_results$y      = factor(df_results$y     ,levels= c("0", "1"), labels=c("No", "Yes"))

head(df_results)
raster::spplot(raster::stack(d, p), main = "Prediction HyMapD")
caret::confusionMatrix(data=as.factor(df_results$y_pred), reference=as.factor(df_results$y), positive="Yes")

p1 <- p
p1[is.na(b[[1]])] <- NA


#########################################
### Binary logistic regression statistics
if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

f<- function(lrn, lrn_model_num){
  cat("Loading learner\n")
  lr_n <- lrn$learners[[lrn_model_num]]
  cat("Loading model\n")
  m_n <- lr_n$model
  test_set_n <- lrn$resampling$test_set(lrn_model_num)
  # train_set_n <- lrn$resampling$train_set(lrn_model_num)
  test_set_n <- lrn$task$data()[test_set_n]
  pred_resp <- predict(m_n, newdata=test_set_n, type="response")
  table(test_set_n$Geothermal, (pred_resp > 0.5)*1, dnn=c("Truth","Predicted"))
  pred <- ROCR::prediction(pred_resp, test_set_n$Geothermal)
  #Get the AUC
  auc_n = unlist(slot(ROCR::performance(pred, "auc"), "y.values"))
  cat(paste0("AUC :", round(auc_n*100, digits=2), "%\n"))
  perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
  perf <- ROCR::performance(pred, measure = "sens", x.measure = "spec")
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  plot(perf) # , colorize=TRUE)
  return (perf)
}

f2<- function(lrn, lrn_model_num){
  cat("Loading learner\n")
  lr_n <- lrn$learners[[lrn_model_num]]
  cat("Loading model\n")
  m_n <- lr_n$model
  test_set_n <- lrn$resampling$test_set(lrn_model_num)
  # train_set_n <- lrn$resampling$train_set(lrn_model_num)
  test_set_n <- lrn$task$data()[test_set_n]
  pred_resp <- 1-predict(m_n, newdata=test_set_n)#, type="response")
  table(test_set_n$Geothermal, (pred_resp > 0.5)*1, dnn=c("Truth","Predicted"))
  pred <- ROCR::prediction(pred_resp, test_set_n$Geothermal)
  #Get the AUC
  auc_n = unlist(slot(ROCR::performance(pred, "auc"), "y.values"))
  cat(paste0("AUC :", round(auc_n*100, digits=2), "%\n"))
  perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
  perf <- ROCR::performance(pred, measure = "sens", x.measure = "spec")
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  plot(perf) # , colorize=TRUE)
  return (perf)
}

pacman::p_load(caret, pROC, ROCR)

perf_glm <- f(rr_lrcv, 4) 
glm_logistic_regression_4 = rr_lrcv$learners[[4]]$model
p <- raster::predict(m, glm_logistic_regression_4)
raster::spplot(raster::stack(m[[1]], p>0), main="LR Prediction>0")

################ 
######### https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
######### https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_mult_logistic_gof_pseudo_r_squared.htm
######### https://www.statology.org/glm-r-squared/
######### McFadden's pseudo R-squared
######### https://scholar.google.com/scholar?cluster=5227738984124726719

McFadden_lr4 = with(summary(glm_logistic_regression_4), 1 - deviance/null.deviance)
cat(paste0("McFadden pseudo R-squared: ", McFadden_lr4, "\n"))
plot(perf) # , colorize=TRUE)
roc_4 <- pROC::roc(predictor=pred_resp, response=test_set_4$Geothermal, na.rm=TRUE, levels=c("No", "Yes"))

install.packages('pscl')
library(pscl)
pscl::pR2(rr_lrcv$learner)['McFadden']
train_set_n <- rr_lrcv$resampling$train_set(4)
glm_logistic_regression_4=glm(Geothermal~., family = "binomial", data = rr_lrcv$task$data()[train_set_n])
pscl::pR2(glm_logistic_regression_4)#['McFadden']
