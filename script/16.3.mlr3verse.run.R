################################################################################
##### setting
library(data.table)
library(openxlsx)
library(mlr3verse); library(dplyr); library(broom); library(factoextra); library(future); 


rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq"; set.seed(100)
dir_dataset = file.path(dir_main, "16.mlr3verse/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_impo = file.path(dir_main, "16.mlr3verse/importance"); if(!dir.exists(dir_impo)) dir.create(dir_impo, recursive = TRUE)

# choose subtype
args <- commandArgs(trailingOnly = TRUE)
choo_cluster <- args[1]

dir_result = file.path(dir_main, paste0("16.mlr3verse/result_", choo_cluster)); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
setwd(dir_result);

cycle = 5  
################################################################################
##### read data
dataset <- readRDS(file.path(dir_dataset, "dataset.rds"))
dataset_data <- dataset
rownames(dataset_data) <- dataset_data[, 1]
dataset_data <- dataset_data[, -1]
colnames(dataset_data) <- gsub("-", "_", colnames(dataset_data))  
colnames(dataset_data) <- gsub("/", ".", colnames(dataset_data)) 
dataset_data$state[dataset_data$state != choo_cluster] <- "Hep_others"

################################################################################
##### target
gene_feature <- readRDS(file.path(dir_dataset, "gene_feature.rds"))
gene_feature <- gene_feature[1:2000]
gene_feature <- gsub("-", "_", gene_feature)  
gene_feature <- gsub("/", ".", gene_feature) 
# print(length(intersect(colnames(dataset_data), gene_feature)))  # 查看共有基因


dataset.task <- dataset_data[, colnames(dataset_data) %in% c("state", gene_feature)]
################################################################################
##### task
task = as_task_classif(dataset.task, target = "state") 

split = partition(task, ratio = 0.8) 
train = task$clone()$filter(split$train)   
test = task$clone()$filter(split$test)   
################################################################################
##### learner
#### pipeline
linear_pipeline = po("scale") %>>%
  po("learner", learner = lrn("classif.ranger", predict_type = "prob"))

glrn = as_learner(linear_pipeline)

################################################################################
##### automatic parameterization
search_space = ps(
  classif.ranger.mtry = p_int(lower = 5, upper = 30),
  classif.ranger.num.trees = p_int(lower = 100, upper = 1500),
  classif.ranger.min.node.size = p_int(lower = 5, upper = 20))

at <- auto_tuner(
  tuner = tnr("random_search"),
  learner = glrn,
  resampling = rsmp("cv", folds = 5), 
  measure = msr("classif.auc"),
  search_space = search_space,
  term_evals = 5
)


at$train(task)  

at$tuning_result   
parameter <- at$tuning_result$learner_param_vals[[1]]
saveRDS(parameter, file.path(dir_result, paste0("parameter_", choo_cluster, ".rds")))


##### hyperparameterization
parameter <- readRDS(file.path(dir_result, paste0("parameter_", choo_cluster, ".rds")))
glrn$param_set$values = parameter  # 更新图学习器超参数


################################################################################
##### Importance
if (T){
  learner_rf = lrn("classif.ranger", importance = "permutation")
  learner_rf$model$learner_params = parameter
  learner_rf$train(task, row_ids = split$train)

  
  importance = as.data.table(learner_rf$importance(), keep.rownames = TRUE)
  colnames(importance) = c("Feature", "Importance")
  write.xlsx(importance, file = file.path(dir_impo, paste0("importance_", choo_cluster, ".xlsx")), row.names = F, col.names = T)
  saveRDS(learner_rf, file.path(dir_result, paste0("importance_", choo_cluster, ".rds")))
  
  
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
  
  png(file.path(dir_impo, paste0("importance_", choo_cluster, ".png")), width = 3000, height = 3000, res = 400, bg = "transparent"); print(p); dev.off()
}





##### Resampling
resampling = rsmp("cv", folds = cycle)
rr = resample(task, glrn, resampling, store_models = TRUE)

rr$score(msr("classif.acc"))        
rr$aggregate(msr("classif.acc"))    
saveRDS(rr, file.path(dir_result, "rr.rds"))
