# Class & Drug sensitivity analysis

# conda activate oncoPredict
# cd /share/home/shh/HCC_snRNA-seq/36.oncoPredict/result
# nohup R --max-ppsize 500000 -f '/share/home/shh/HCC_snRNA-seq/script/36.2.oncoPredict.run.R' > 36.oncoPredict.log 2>&1 &
################################################################################
##### setting
library(oncoPredict)

rm(list = ls());
dir_main = "/share/home/shh/HCC_snRNA-seq";  set.seed(100)
dir_in = file.path(dir_main, "114.class")
dir_dataset = file.path(dir_main, "36.oncoPredict/dataset"); if(!dir.exists(dir_dataset)) dir.create(dir_dataset, recursive = TRUE)
dir_result = file.path(dir_main, "36.oncoPredict/result"); if(!dir.exists(dir_result)) dir.create(dir_result, recursive = TRUE)
# writeLines(capture.output(sessionInfo()), file.path(dir_dataset, "session.txt"))
setwd(dir_result);
################################################################################
##### data 
CTRP2_Expr = readRDS(file.path(dir_dataset, "Training Data/CTRP2_Expr (TPM, not log transformed).rds"))
CTRP2_Expr[1:4,1:5]

CTRP2_Res = readRDS(file.path(dir_dataset, "Training Data/CTRP2_Res.rds"))
CTRP2_Res[1:4,1:5]

Trai_Expr <- CTRP2_Expr
Trai_Res <- CTRP2_Res


Exp <- readRDS(file.path(dir_in, "dataset/df_ComBat_all.rds"))
Exp <- as.matrix(Exp)
Exp[1:4,1:5]


calcPhenotype(trainingExprData = Trai_Expr,
              trainingPtype = Trai_Res,
              testExprData = Exp,
              batchCorrect = 'eb', 
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')


