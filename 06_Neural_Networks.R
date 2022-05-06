if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(sf, terra, tidyverse, neuralnet, caret)

input_stack_file = "d:/geoai_som/brady_som_output.gri"
s <- terra::rast(input_stack_file)
cat("Stack loaded:")
print(s)
s <- s[[c("Geothermal", "Minerals", "Temperature", "Faults")]]
stack_df <- as.data.frame(s, xy = FALSE)  # we don't need the (x, y) coordinates
samplesize = 0.7*nrow(stack_df)
set.seed(42)
index = sample(seq_len(nrow(stack_df)), size = samplesize)

train_data = stack_df[index, ]
test_data = stack_df[-index, ]
# Not needed. Data is already scaled
# max_sf = apply(stack_df, 2, max)
# min_sf = apply(stack_df, 2, min)
# avg_sf = apply(stack_df, 2, mean)
# stack_df = as.data.frame(scale(stack_df, center = FALSE, scale = max_sf-min_sf))
set.seed(42)
ann_1 = neuralnet::neuralnet(Geothermal ~ Minerals+Temperature+Faults,
                             train_data,
                             threshold=3.0,
                             hidden = 1, # rep = 5,
                             algorithm = 'rprop+', # 'backprop', # or 'rprop+'
                             learningrate = 0.001,
                             lifesign = 'full',
                             lifesign.step = 100,
                             act.fct = 'logistic',
                             stepmax = 10000,
                             linear.output = F)
print(ann_1$result.matrix)
plot(ann_1)


predict_testNN = neuralnet::compute(ann_1, test_data[,c("Minerals", "Temperature", "Faults")])
binary_testNN = round(predict_testNN$net.result, digits=0)
# View(predict_testNN)
# View(test_data$Geothermal)
results = data.frame(actual=test_data$Geothermal, pred=predict_testNN$net.result)
confusion_1 <- caret::confusionMatrix(reference=test_data$Geothermal,
                                      data=predict_testNN$net.result)
                                        
roundedresults=sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
confusion_ann_1 <- as.data.frame(unclass(table(roundedresultsdf$actual,roundedresultsdf$pred)))
accuracy_ann_1 <- (confusion_ann_1[1,1]+confusion_ann_1[2,2])/nrow(roundedresultsdf)
print(confusion_ann_1)
print(accuracy_ann_1)
confusion_1 <- caret::confusionMatrix(reference=
                                        factor(if_else(roundedresultsdf$actual==1, "true", "false"),
                                               levels = c("true", "false")), 
                                      data = factor(if_else(roundedresultsdf$pred==1, "true", "false"),
                                                    levels = c("true", "false")))
print(confusion_1)


ann_2 = neuralnet::neuralnet(Geothermal ~ Minerals+Temperature+Faults,
                             train_data,
                             threshold=4.0, # 0.01,
                             hidden = 4,
                             # rep=5,
                             algorithm = 'rprop+', # 'backprop', # or 'rprop+'
                             learningrate = 0.01,
                             lifesign = 'full',
                             lifesign.step = 100,
                             act.fct = 'logistic',
                             stepmax = 10000,
                             linear.output = F)
plot(ann_2)
