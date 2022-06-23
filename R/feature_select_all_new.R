

# 根据 feature 找 gene 表达数据 --------------------------------------------------
rm(list = ls())

library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用



setwd('D:\\E\\博士\\R_程序\\EFS\\Results')
rank = read.csv(file = "rank.csv", header = T, sep = ",")
gene <- as.matrix(rank[rank[,2] >= 0.4,1])


setwd('D:\\E\\博士\\R_程序\\EFS\\Data')
myfile = list.files("GEO")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./GEO/", myfile, sep = "")     #用paste命令构建路径变量dir
n = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
for (i in 1:n) {
  # i <- 1
  Data1 = read.table(file = dir[i], header = T, check.names = FALSE)
  dim(Data1)    # 21836   185
  # View(Data1[,1:10])
  
  colnames(gene) <- c('gene')
  # Data1 <- t(Data)
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  # View(Data2[,1:10])
  
  genedata <- merge(gene, Data2, by = "gene")
  # View(genedata[,1:10])
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[1,], genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  # View(genedata2)
  
  # i <- 4
  name <- dir[i]
  name1 <- str_split_fixed(name, "./", 2)
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2)
  name4 <- name3[1]
  name5 <- str_c("Feature13_", name4)
  
  
  path <- paste("./FeatureGEO/", paste(name5, ".txt", sep = ""), sep = "")
  write.table(genedata2, path, quote = F, sep = "\t")
  
  
  data1 = read.table(
    "D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA\\TCGA_pro_outcome_TN_log_scale.txt",
    header = T,
    check.names = FALSE
  )
  dim(data1)    # 824 224
  gene <- as.matrix(genedata[, 1])
  
  colnames(gene) <- c('gene')
  # Data1 <- t(Data)
  data2 <- cbind(rownames(data1), data1)
  colnames(data2) <- c('gene', colnames(data1))
  # View(data2[,1:10])
  
  genedata <- merge(gene, data2, by = "gene")
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(data1[1,], genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  
  name5 <- str_c("TCGA_", name4)
  
  
  path <- paste("./FeatureTCGA/", paste(name5,".txt", sep=""), sep="")
  write.table(genedata2, path, quote = F, sep="\t")
  
}              





# classo performance  -----------------------------------------------------

## clear
rm(list = ls())

## package
library(pROC)
library(stringr)
myname <- function(x){
  # x <- 1
  name <- dir[x]
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[1], "[_]", 2);
  name4 <- str_split_fixed(name3[2], "[.]", 2);
  name5 <- str_split_fixed(name4[1], "[_]", 2);
  name6 <- str_c("feature_",name5[1])
  return(name6)
}


## 分类指标的函数
classindex <- function(A_test, Data, thro){
  predict = ifelse(A_test[,1] >= thro, 1, 0)
  predict_value = predict
  true_value = A_test[,2]
  error = predict_value-true_value
  
  data <- t(Data)
  # 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
  accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
  precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
  recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
  # P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
  F_measure = 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
  specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
  table = table(predict_value, true_value) 
  result <- list(accuracy, precision, recall, specificity, F_measure,table)
  names(result) <- c("accuracy", "precision", "recall", "specificity", "F_measure", "table")
  return(result)
}

## input data
setwd("D:\\E\\博士\\R_程序\\EFS\\Data")


myfile = list.files("FeatureGEO")    
dir = paste("FeatureGEO/", myfile, sep="")   
n = length(dir) 


auc <- c()
for (i in 1:n) {
  # i <- 14
  myfile = list.files("FeatureGEO")    
  dir = paste("FeatureGEO/", myfile, sep="") 
  Data  = read.table(file = dir[i], header = T, check.names = FALSE)
  
  myfile = list.files("FeatureTCGA")    
  dir = paste("FeatureTCGA/", myfile, sep="") 
  Data2 = read.table(file = dir[i], header = T, check.names = FALSE)
  
  x.train <- data.frame(t(Data2)[,-1])
  y.train <- t(Data2)[,1]
  x.test <- data.frame(t(Data)[,-1])
  y.test <- t(Data)[,1]
  
  # model <- svm(x.train, y.train, kernel = 'linear', scale = FALSE)   # linear  radial
  # summary(model)
  
  model <- glm(y.train~., data = x.train, family = binomial)
  # glmfit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
  summary(model)
  
  ## prediction
  p_test <- predict(model, x.test, type = "response")
  p_test = as.matrix(p_test)
  A_test <- data.frame(p_test, y.test)
  names(A_test)<- c("p", "outcome")
  
  p <- A_test[,1]
  p_glm <- cbind(log(p/(1-p)), A_test[,2])
  colnames(p_glm) <- c('y.test', 'Lable')
  
  path <- paste("./FeaturePRS/",paste(myname(i),".csv"))
  write.csv(p_glm, path, row.names = F)

  
  path <- paste("./FeatureAUC/",paste(myname(i),".csv"))
  write.csv(A_test, path, row.names = F)
  index <- classindex(A_test, t(x.test), 0.500)
  picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
  auc <- cbind(auc, picauc$auc)
  print(picauc$auc)
  # print(index)
  print(i)
}



