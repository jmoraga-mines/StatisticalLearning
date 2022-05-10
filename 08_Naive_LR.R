if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC)

# "d:/ThesisLayers/AI_Stacks/FullStacks/desert_full_stack"
# "d:/ThesisLayers/AI_Stacks/FullStacks/brady_full_stack"
# "d:/ThesisLayers/AI_Stacks/MineralStacks/desert_minerals_stack"
# "d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack"
m <- raster::stack("d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack")

new_ai_data_df <- as.data.frame(m, xy=F)
new_ai_data_df <- na.omit(new_ai_data_df)
# new_ai_data_df$Geothermal <- as.factor(new_ai_data_df$Geothermal)

new_ai_data_df$Geothermal <- factor(ifelse(new_ai_data_df$Geothermal>=0.5, "Yes", "No"),
                                    levels = c("No", "Yes"))

true_data <- new_ai_data_df[new_ai_data_df$Geothermal=="Yes", ]
false_data <- new_ai_data_df[new_ai_data_df$Geothermal=="No", ]
true_data <- unique(true_data)
false_data <- unique(false_data)
set.seed(42)
true_idx <- sample(seq_len(nrow(true_data)), size = 20000)
false_idx <- sample(seq_len(nrow(false_data)), size = 20000)
train_data <- rbind(true_data[true_idx, ], false_data[false_idx, ])
test_data <- rbind(true_data[-true_idx, ], false_data[-false_idx, ])
nrow(train_data)
nrow(test_data)

set.seed(123)
lr <- glm(Geothermal ~ .,  
          family = binomial, data = train_data)
set.seed(123)
lr <- glm(Geothermal ~ .,  #  -Epsomite-Minerals-Deformation-Chalcedony
          family = binomial, data = new_ai_data_df) # Or if factors use binomial

print(car::vif(lr))
geo_predict <- predict(lr, newdata=test_data, type="response") # "class" or "response"
geo_predict_factor <- geo_predict
geo_predict_factor <- factor(ifelse(geo_predict>=0.5, "Yes", "No"), levels = c("No", "Yes"))
c <- confusionMatrix(data = geo_predict_factor, 
                     reference = test_data$Geothermal,
                     positive = "Yes")
print(c)
pROC::auc(test_data$Geothermal, geo_predict)
geo_roc <- pROC::roc(test_data$Geothermal, geo_predict)

r <- raster::predict(m, lr)
r <- r>0.5
names(r) <- "LR"
raster::spplot(raster::stack(m[["Geothermal"]], r), main="Prediction")
