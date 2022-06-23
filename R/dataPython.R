
rm(list = ls())


## 原始数据
setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA')

Data <- read.table("TCGA_pro_outcome_TN_log_train.txt", header = T, check.names = FALSE)
data <- read.table("TCGA_pro_outcome_TN_log_trainP.txt", header = T, check.names = FALSE)

dim(data)    # 158 824
# View(data[,1:10])

## Python data
setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')


# varicance ---------------------------------------------------------------

var <- read.table("variance.txt", header = T, check.names = FALSE)
rownames(var) <- rownames(data)

data[,2]
var[,2]

lab <- c()
for (i in 1:dim(data)[2]) {
  for (j in 1:dim(var)[2]) {
    if (sum(data[,i] - var[,j]) == 0)
      lab <- rbind(lab, i)
  }
}

## 没选出的为0
imp = cbind(1:dim(data)[2],  rep(0, dim(data)[2])) 
imp[lab ,2] <- 1

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankvar = imp[,1]
impvar = imp[,2]


if(length(imp1)>norm){
  impvar = impvar[1:norm]
}
Variance <-  impvar[order(rankvar, decreasing = F)]



# chi2 --------------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')
chi2 <- read.table("chi2.txt", header = T, check.names = FALSE)


imp = cbind(1:dim(data)[2], chi2/max(chi2)) 

## 设置阈值，定义重要性
imp[which(imp[,2] <= 0.3),2] <- 0

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankchi2 = imp[,1]
impchi2 = imp[,2]


if(length(impchi2)>norm){
  impchi2 = impchi2[1:norm]
}
chi2 <-  impchi2[order(rankchi2, decreasing = F)]




# mi ----------------------------------------------------------------------

setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')
mi <- read.table("mi.txt", header = T, check.names = FALSE)

imp = cbind(1:dim(data)[2], mi/max(mi)) 

## 设置阈值，定义重要性
imp[which(imp[,2] <= 0.5),2] <- 0

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankmi = imp[,1]
impmi = imp[,2]


if(length(impmi)>norm){
  impmi = impmi[1:norm]
}
mi <-  impmi[order(rankmi, decreasing = F)]



# Relief ------------------------------------------------------------------

setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')
relief <- read.table("relief.txt", header = T, check.names = FALSE)

imp = cbind(1:dim(data)[2], abs(relief)/max(abs(relief))) 
## 设置阈值，定义重要性
imp[which(imp[,2] <= 0.1),2] <- 0

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankrelief = imp[,1]
imprelief = imp[,2]


if(length(imprelief)>norm){
  imprelief = imprelief[1:norm]
}
relief <-  imprelief[order(rankrelief, decreasing = F)]


# DT ------------------------------------------------------------------

setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')
dt <- read.table("dt.txt", header = T, check.names = FALSE)

imp = cbind(1:dim(data)[2], dt[,1]/max(dt[,1])) 
## 设置阈值，定义重要性
imp[which(imp[,2] <= 0.5),2] <- 0

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankdt = imp[,1]
impdt = imp[,2]


if(length(impdt)>norm){
  impdt = impdt[1:norm]
}
dt <-  impdt[order(rankdt, decreasing = F)]




# mRMR --------------------------------------------------------------------

rm(list = ls())

# install.packages("mRMRe")
library(mRMRe)
library(caret)

seed <- 123    # svmfs 函数必须的

setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA') 

## 总数据
data_all <- read.table("TCGA_pro_outcome_TN_log_train.txt", header = T, check.names = FALSE, sep="\t")
View(data_all[,1:10])
dim(data_all)    # 158 824
data <- data.frame(data_all)
View(data[,1:10])


data[,1] <- factor(data[,1])
data[,1] <- ordered(data[,1], levels = c("0", "1"))


x <- data[,-1]    # sample*gene
View(x[,1:10])
y <-  as.numeric(as.factor(data[,1]))-1 


data.mrmre.train <- mRMR.data(data=data[,-1], strata = data[,1])
## classical mRMR feature selection
# res.fs.mrmr <- mRMR.classic(data=data.mrmre.train, target_indices=1, feature_count=25)
## 与上等价：solution_count不是传统mRMR的参数，所以如果设置为1，它仅执行一个经典的mRMR
res.fs.mrmr <- mRMR.ensemble(data.mrmre.train, target_indices=1, feature_count=200, solution_count = 1) # , method="exhaustive"
list(apply(solutions(res.fs.mrmr)[[1]], 2, function(x, y) 
{ return(y[x]) }, y = featureNames(data.mrmre.train)))
# myfs <- c(myfs, list(apply(solutions(res.fs.mrmr)[[1]], 2, function(x, y) 
# { return(y[x]) }, y = featureNames(data.mrmre.train))))
selected.features.mrmre <- mRMRe::solutions(res.fs.mrmr)
class(selected.features.mrmre)



## 没选出的为0
imp = cbind(1:dim(data)[2],  rep(0, dim(data)[2]))
lab <- selected.features.mrmre$`1`

imp[lab, 2] <- 1

imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rankmRMR = imp[,1]
impmRMR  = imp[,2]


if(length(impmRMR)>norm){
  impmRMR = impmRMR[1:norm]
}
mRMR <-  impmRMR[order(rankmRMR, decreasing = F)]








