# conda activate class && R
################################################################################
##### setting
library(dplyr)
library(data.table)
library(openxlsx)
library(mlr3verse); library(dplyr); library(broom); library(factoextra); library(future); 

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_in = file.path(dir_main, "16.mlr3verse/importance")
dir_dataset = file.path(dir_main, "114.class/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_mlr3 = file.path(dir_main, "114.class/dataset/mlr3verse"); if(!dir.exists(dir_mlr3)) dir.create(dir_mlr3, recursive = TRUE)
dir_result = file.path(dir_main, "114.class/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result)

cycle = 5
################################################################################
##### read data
df_Class <- readRDS(file.path(dir_main, '114.class/result', paste0('df_Class.rds')))
df_exp <- readRDS(file.path(dir_dataset, paste0("df_ComBat.rds")))
gene_feature <- readRDS(file.path(dir_dataset, paste0("gene_feature.rds"))) 
df_exp <- df_exp[gene_feature,]
df_trai <- df_exp[,colnames(df_exp) %in% df_Class$Sample]
df_trai <- df_trai[,match(df_Class$Sample, colnames(df_trai))]

# state
state <- df_Class$Class %>% data.frame()  # 细胞类型
rownames(state) <- df_Class$Sample
colnames(state) <- 'state'
# table(state$state)

data_tmp <- t(df_trai) %>% data.frame()
data_trai <- merge(data_tmp, state, by = "row.names")
rownames(data_trai) <- data_trai[, 1]
data_trai <- data_trai[, -1]


dataset.task <- data_trai
# write.csv(dataset.task, file = "dataset_task.csv", row.names = FALSE)

################################################################################
##### task
task <- as_task_classif(dataset.task, target = "state") 
split = partition(task, ratio = 0.8) 
train = task$clone()$filter(split$train)    
test = task$clone()$filter(split$test)     

################################################################################
##### learner
linear_pipeline = po("scale") %>>%
  po("learner", learner = lrn("classif.ranger", predict_type = "prob"))

glrn = as_learner(linear_pipeline)

################################################################################
##### auto_tuner
if(T){
  at <- auto_tuner(
    tuner = tnr("random_search"),
    learner = glrn,
    resampling = rsmp("cv", folds = 3), 
    measure = msr("classif.acc"),
    term_evals = 20 
  )
  
  at$train(task, row_ids = split$train) 
  
  at$tuning_result           
  parameter <- at$tuning_result$learner_param_vals[[1]]
  saveRDS(parameter, file.path(dir_mlr3, paste0("parameter_trai.rds")))
}

parameter <- readRDS(file.path(dir_mlr3, paste0("parameter_trai.rds")))
glrn$param_set$values = parameter  # 更新图学习器超参数


################################################################################
##### importance
if (T){
  learner_rf = lrn("classif.ranger", importance = "permutation")
  learner_rf$model$learner_params = parameter
  learner_rf$train(task, row_ids = split$train)

  
  importance = as.data.table(learner_rf$importance(), keep.rownames = TRUE)
  colnames(importance) = c("Feature", "Importance")
  write.xlsx(importance, file = file.path(dir_mlr3, paste0("importance_trai.xlsx")), row.names = F, col.names = T)
  saveRDS(learner_rf, file.path(dir_mlr3, paste0("importance_trai.rds")))
  
  
  p <- ggplot(data = importance[1:20], aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col() +
    coord_flip() +
    xlab("") +
    ylab("coefficient") +
    ggtitle("Importance (Top20)") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          plot.title = element_text(size = 30, hjust = 0.5))
  
  png(file.path(dir_result, paste0("114.importance_trai.png")), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()
}


################################################################################
##### train
if (T){
  glrn$train(task, row_ids = split$train)
  glrn$model$classif.ranger   
  saveRDS(glrn, file.path(dir_mlr3, paste0("glrn_trai.rds")))

  predictions = glrn$predict(task, row_ids = split$test)
  predictions$score(msr("classif.acc"))     
  saveRDS(predictions, file.path(dir_mlr3, paste0("predictions_trai.rds")))
}


################################################################################
##### predict
df_vali <- df_exp[,!colnames(df_exp) %in% df_Class$Sample]
new_data <- t(df_vali) %>% data.frame()

tmp_vali <- glrn$predict_newdata(new_data)
pred_vali <- data.frame(Sample = rownames(new_data) , Class = tmp_vali$response)
table(pred_vali$Class)
table(data_trai$state)
saveRDS(pred_vali, 'pred_vali.rds')

