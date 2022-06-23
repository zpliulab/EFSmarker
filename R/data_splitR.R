
rm(list = ls())

library(caret)
library(dplyr)

# 路径 ----------------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA')
Data <- read.table("TCGA_pro_outcome_TN_log_UN824.txt",header=T,sep='\t', check.names = F)
View(Data[,1:10])
dim(Data)   # 824   224


##  scale
all_data_scale <- rbind(Data[1,], t(scale(t(Data[-1,]))))
dim(all_data_scale)    # 824   224
View(all_data_scale[,1:10])
# write.table(all_data_scale, file = "TCGA_pro_outcome_TN_log_scale.txt",quote = F, sep = "\t")



# 数据 ----------------------------------------------------------------------
# x <- read.table("TCGA_pro_outcome_TN_log_scale.txt", header = T, check.names = FALSE)
x <- Data
data <- data.frame(t(x))
View(data[,1:10])

# 均值和方差
i <- 12
mean(data[,i])
var(data[,i])

# 数据 —— 训练集+测试集 -----------------------------------------------------------

## 第一次试验的结果
set.seed(123*1)
training.samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- as.matrix(data[training.samples, ])
test.data <- as.matrix(data[-training.samples, ])  


write.table(train.data, file = "TCGA_pro_outcome_TN_log_train.txt",quote = F, sep = "\t")
write.table(test.data, file = "TCGA_pro_outcome_TN_log_test.txt",quote = F, sep = "\t")

