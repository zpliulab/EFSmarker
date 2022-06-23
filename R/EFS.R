

rm(list = ls())

## 
library(EFS)
library(glmnet)


setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA')


datatrain <- read.table("TCGA_pro_outcome_TN_log_train.txt", header = T, check.names = FALSE)
datatest <- read.table("TCGA_pro_outcome_TN_log_test.txt", header = T, check.names = FALSE)
dim(datatrain)    # 158 824
# View(datatrain[,1:10])


### initial

data <- datatrain
classnumber <- 1
NA_threshold = 0.2 
cor_threshold = 0.7
runs = 2
# selection = c(TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, FALSE, FALSE)
# selection = c(TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE, TRUE)


start.time <- Sys.time()
classname = colnames(data)[classnumber]


# L??schen von Parametern mit zu vielen NA
NrNA= c()
for(i in 1:length(data[1,])){
  NrNA[i] = length(which(is.na(data[,i])))
}
NmNA = which(NrNA/length(data[,1])>NA_threshold)
if(length(NmNA) != 0) data=data[,-NmNA]

data = na.omit(data, row.names=F)

# oder nur Nullen oder Varianz gleich Null
data=data[,which(colSums(data,na.rm=F)!=0 &
                   apply(data, 2, var)!=0)]

klasse = data[,classname]
# teste auf bin??re Klasse
is.binary <- function(v) {
  x <- unique(v)
  length(x) == 2L && all(x==1 | x==0)
}
if(!is.binary(klasse)){
  stop("Class not binary with classlabels 0 and 1")
}

clnr= which(colnames(data)==classname)
data = data[,-clnr]
norm=length(data[1,])




# median filter -----------------------------------------------------------
print('Start Median')
positives = which(klasse==1)
negatives = which(klasse==0)

data.pos = data[positives,]
data.neg = data[negatives,]


f <- function(x,y){
  test <- wilcox.test(x,y)
  data.frame(pval = test$p.value)
}
wtest=sapply(colnames(data), function(x) f(data.pos[,x],
                                           data.neg[,x]))
## p value
pval <- unlist(wtest, use.names=FALSE)
imp = cbind(1:length(pval),  1-pval+min(pval) )
# colnames(data)[1:10]
# View(imp[,2])
## 将重要性 <1 的设为0
imp[which(imp[,2] < 1), 2] <- 0


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[,which(imp[,2] < 1)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[,which(imp[,2] < 1)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
library(pROC)
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucmedian <- picauc$auc


# coef <- as.matrix(model$coefficients)
# coef00 <- coef[-1,1]
# lab <- which(coef00 != 'NA')
# coef00[-lab] <- 0
# coef11 <- rbind(coef[1,1], as.matrix(coef00))


imp = imp[order(imp[,2], decreasing = T),]
# order(imp[,2])   # incrsing
## 重要性从高到低排序
rank1 = imp[,1]
imp1 = imp[,2]*aucmedian 

if(length(imp1)>norm){
  imp1 = imp1[1:norm]
}
Median <-  imp1[order(rank1, decreasing = F)]
length(which(Median!=0))


medianwei <- cbind(colnames(data), Median)
# write.csv(medianwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/medianfea.csv", row.names = F)
medianfea <- colnames(data)[which(Median!=0)]
# write.csv(medianfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/medianfea.csv", row.names = F)



# Pearson correlation filter ----------------------------------------------
cor.filter <- function(data, threshold, method){
  
  #catch missing or bad input
  
  if(!is.numeric(threshold) | threshold > 1 |
     threshold < 0)
    stop("invalid argument:
         threshold is required to be in [0,1]")
  
  #analysis
  # data.tmp = as.matrix(cbind(klasse, data))
  data.tmp = data
  features = c()
  
  klasse = data.tmp[,1]
  data.tmp = data.tmp[,-1]
  colnames(data.tmp) = seq(1,length(data.tmp[1,]))
  
  vars = apply(data.tmp, 2, var)
  dels = which(vars == 0)
  if(length(dels) > 0){
    data.tmp = data.tmp[,-dels]
  }
  
  names = colnames(data.tmp)
  
  for(k in 1:norm){
    # k <- 1
    #h??chste Korrelation mit Class label
    # c.matrix = cor(cbind(klasse, data.tmp), method="pearson")
    c.matrix = cor(cbind(klasse, data.tmp), method=method)
    # View(c.matrix[,1:10])
    
    # Betrag der ersten Spalte ohne erste Zeile
    #korrelation mit class
    cor.class = abs(as.vector(c.matrix[1,-1]))
    
    #transponieren und als mit
    #names = seq(1, length(data.temp[1,])) versehen
    
    imp = t(as.matrix(cor.class))
    best = names[which.max(as.vector(imp[1,]))]
    
    features = c(features, best)
    
    #########
    
    c.matrix = c.matrix[-1,-1]
    tmp = abs(c.matrix[best,])
    
    # w = which(tmp > cor_threshold)
    w = which(tmp > threshold)
    
    
    data.tmp = data.tmp[,-w]
    
    names = names[-w]
    
    if(is.null(dim(data.tmp))){
      features = as.numeric(features)
      return(features)
    }
    colnames(data.tmp) = names
  }
  return(features)
}

print('Start Pearson')
rank2 = cor.filter(as.matrix(cbind(klasse, data)),
                   cor_threshold, "pearson")
if(length(rank2)>norm){
  rank2 = rank2[1:norm]
}


c.matrix = cor(cbind(klasse, data), method="pearson")
cor.class = abs(as.vector(c.matrix[1,-1]))
imp = t(as.matrix(cor.class))
imp=as.vector(imp)
imp= imp[rank2]
imp2 = imp/max(imp)

#dependent variables omitted
#refill vector to length = norm
if(length(imp2)>norm){
  imp2 = imp2[1:norm]
}

if(length(rank2)<norm){
  rest=setdiff(c(1:norm),rank2)
  k=length(rest)
  rank2=c(rank2,rest)
  imp2=c(imp2,rep(0,k))
}

P_cor <- imp2[order(rank2, decreasing = F)]
# View(P_cor)


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(P_cor == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(P_cor == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucPcor <- picauc$auc

P_cor <- imp2[order(rank2, decreasing = F)]*aucPcor
length(which(P_cor!=0))


P_corwei <- cbind(colnames(data), P_cor)
# write.csv(P_corwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/P_corwei.csv", row.names = F)
P_corfea <- colnames(data)[which(P_cor!=0)]
# write.csv(P_corfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/P_corfea.csv", row.names = F)



# Pearson correlation filter ----------------------------------------------
print('Start Spearman')
rank3 = cor.filter(as.matrix(cbind(klasse, data)),
                   cor_threshold, "spearman")
if(length(rank3)>norm){
  rank3 = rank3[1:norm]
}

c.matrix = cor(cbind(klasse, data), method="spearman")
cor.class = abs(as.vector(c.matrix[1,-1]))
imp = t(as.matrix(cor.class))
imp=as.vector(imp)
imp= imp[rank3]
imp3 = imp/max(imp)

#dependent variables omitted
#refill vector to length = norm
if(length(imp3)>norm){
  imp3 = imp3[1:norm]
}

if(length(rank3)<norm){
  res=setdiff(c(1:norm),rank3)
  l=length(res)
  rank3=c(rank3,res)
  imp3=c(imp3,rep(0,l))
}

S_cor <- imp3[order(rank3, decreasing = F)]



## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(S_cor == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(S_cor == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucScor <- picauc$auc

S_cor <- imp3[order(rank3, decreasing = F)]*aucScor
length(which(S_cor!=0))


S_corwei <- cbind(colnames(data), S_cor)
# write.csv(S_corwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/S_corwei.csv", row.names = F)
S_corfea <- colnames(data)[which(S_cor!=0)]
# write.csv(S_corfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/S_corfea.csv", row.names = F)





# Ridge regression -----------------------------------------------------

print('Start Ridge')
library(glmnet)
imp_ridge <- function(dataf){
  # Z-Transformation:
  dataf <- as.matrix(data)
  # View(dataf[,1:10])
  # class(dataf)
  for (i in 1:length(dataf[1,])){
    ## N(0,1)
    dataf[,i]=(dataf[,i]-mean(dataf[,i]))/var(dataf[,i])
  }
  # klassef = (klasse-mean(klasse))/var(klasse)
  
  ## cv
  cv.ridge = cv.glmnet(dataf, klasse, alpha = 0, family = "binomial", nfolds = 10, type.measure = "class")
  # plot(cv.ridge)
  ## fit
  ridge.model <- glmnet(dataf, klasse, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.1se)
  # plot(ridge.model)
  corri = abs(as.vector(coef(ridge.model)))
  
  # model.imp = glm(as.factor(klassef)~., data = dataf,
  #                 family = binomial(link = "logit"), maxit=100)
  # corri = abs(as.vector(model.imp$coefficients))
  
  corri <- corri[-1]
  corri <- na.omit(corri)
  # View(corri)
  imp <- c()
  imp <- cbind(1:length(corri),  (corri*1)/max(corri))
  # View(imp)
  
  ## 设置阈值，定义重要性
  imp[which(imp[,2] <= 0.5),2] <- 0
    
  
  imp <- imp[order(imp[,2], decreasing = TRUE),]
  rank4 <- imp[,1]
  imp4 <- imp[,2]
  
  #NA omitted - refill vector to length = norm
  if(length(rank4)<norm){
    k <- norm-length(rank4)
    rest <- setdiff(c(1:norm),rank4)
    rank4 <- c(rank4,rest)
    imp4 <- c(imp4,rep(0,k))}
  
  if(length(imp4)>norm){
    imp4 <- imp4[1:norm]
  }
  impra4 <- cbind(imp4,rank4)
  return(impra4)
}
impra4 <- imp_ridge(data)

Ridge <- impra4[,1][order(impra4[,2], decreasing = F)]



## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(Ridge == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(Ridge == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucRidge <- picauc$auc

Ridge <- impra4[,1][order(impra4[,2], decreasing = F)]*aucRidge
# View(Ridge)
length(which(Ridge!=0))


Ridgewei <- cbind(colnames(data), Ridge)
# write.csv(Ridgewei, file = "D:/E/博士/R_程序/EFS/Results/Weight/Ridgewei.csv", row.names = F)
Ridgefea <- colnames(data)[which(Ridge!=0)]
# write.csv(Ridgefea, file = "D:/E/博士/R_程序/EFS/Results/Feature/Ridgefea.csv", row.names = F)


# RandomForest - Breiman --------------------------------------------------
print('Start RF')
imp_rf <- function(runs,classnumber){
  
  data.frame= cbind(data,klasse)
  
  imp_accuracy = c()
  imp_Gini = c()
  
  for(i in 1:runs){
    print(i)
    start.run <- Sys.time()
    library(randomForest)
    rf =randomForest(as.factor(klasse)~.,
                     data=data.frame,importance=TRUE,
                     replace=FALSE,ntree=1000)
    imp_accuracy = cbind(imp_accuracy,rf$importance[,3])
    imp_Gini = cbind(imp_Gini,rf$importance[,4])
    end.run <- Sys.time()
    diff <- end.run - start.run
    print(diff)
  }
  
  # Mitteln der runs
  imp_Gini_mean = c()
  imp_Gini_mean = rowMeans(imp_Gini)
  imp_accuracy_mean = c()
  imp_accuracy_mean = rowMeans(imp_accuracy)
  
  #Ranking
  imp= rbind(abs(imp_accuracy_mean),abs(imp_Gini_mean))
  rank_RF1=cbind(1:length(imp[1,]),as.vector(imp[1,]))
  rank_RF1=rank_RF1[order(rank_RF1[,2], decreasing = T),]
  imp_RF1 =rank_RF1[,2]
  rank_RF1=rank_RF1[,1]
  rank_RF2=cbind(1:length(imp[1,]),as.vector(imp[2,]))
  rank_RF2=rank_RF2[order(rank_RF2[,2], decreasing = T),]
  imp_RF2 =rank_RF2[,2]
  rank_RF2=rank_RF2[,1]
  #Normieren
  imp_RF1 = imp_RF1/max(imp_RF1)
  imp_RF2 = imp_RF2/max(imp_RF2)
  imp_RF=cbind(imp_RF1,rank_RF1,imp_RF2,rank_RF2)
  
  return(imp_RF)
}
imp_RF=imp_rf(runs, classnumber)

# RF_ACC -----------------------------------------------------------------

RF_Acc = imp_RF[,1][order(imp_RF[,2], decreasing = F)]


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(RF_Acc == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(RF_Acc == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucRF_Acc <- picauc$auc

RF_Acc <- imp_RF[,1][order(imp_RF[,2], decreasing = F)]*aucRF_Acc
length(which(RF_Acc!=0))


RF_Accwei <- cbind(colnames(data), RF_Acc)
# write.csv(RF_Accwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/RF_Accwei.csv", row.names = F)
RF_Accfea <- colnames(data)[which(RF_Acc!=0)]
# write.csv(RF_Accfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/RF_Accfea.csv", row.names = F)


# RF_Gini -----------------------------------------------------------------

RF_Gini = imp_RF[,3][order(imp_RF[,4], decreasing = F)]

## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(RF_Gini == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(RF_Gini == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucRF_Gini <- picauc$auc

RF_Gini <- imp_RF[,3][order(imp_RF[,4], decreasing = F)]*aucRF_Gini
length(which(RF_Gini!=0))


RF_Giniwei <- cbind(colnames(data), RF_Gini)
# write.csv(RF_Giniwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/RF_Giniwei.csv", row.names = F)
RF_Ginifea <- colnames(data)[which(RF_Gini!=0)]
# write.csv(RF_Ginifea, file = "D:/E/博士/R_程序/EFS/Results/Feature/RF_Ginifea.csv", row.names = F)


# varicance ---------------------------------------------------------------
## 原始数据
setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA')
data <- read.table("TCGA_pro_outcome_TN_log_trainP.txt", header = T, check.names = FALSE)
dim(data)    # 158 824
# View(data[,1:10])

## Python data
setwd('D:\\E\\博士\\R_程序\\EFS\\DataPython')

var <- read.table("variance.txt", header = T, check.names = FALSE)
rownames(var) <- rownames(data)

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


if(length(impvar)>norm){
  impvar = impvar[1:norm]
}
Variance <-  impvar[order(rankvar, decreasing = F)]
# View(Variance)

## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(Variance == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(Variance == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucVariance <- picauc$auc

Variance <- impvar[order(rankvar, decreasing = F)]*aucVariance
length(which(Variance!=0))


Variancewei <- cbind(colnames(data), Variance)
# write.csv(Variancewei, file = "D:/E/博士/R_程序/EFS/Results/Weight/Variancewei.csv", row.names = F)
Variancefea <- colnames(data)[which(Variance!=0)]
# write.csv(Variancefea, file = "D:/E/博士/R_程序/EFS/Results/Feature/Variancefea.csv", row.names = F)




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

Chi2 <-  impchi2[order(rankchi2, decreasing = F)]


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(Chi2 == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(Chi2 == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucChi2 <- picauc$auc

Chi2 <- impchi2[order(rankchi2, decreasing = F)]*aucChi2
length(which(Chi2!=0))


Chi2wei <- cbind(colnames(data), Chi2)
# write.csv(Chi2wei, file = "D:/E/博士/R_程序/EFS/Results/Weight/Chi2wei.csv", row.names = F)
Chi2fea <- colnames(data)[which(Chi2!=0)]
# write.csv(Chi2fea, file = "D:/E/博士/R_程序/EFS/Results/Feature/Chi2fea.csv", row.names = F)



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

MI <-  impmi[order(rankmi, decreasing = F)]


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(MI == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(MI == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucMI <- picauc$auc

MI <- impmi[order(rankmi, decreasing = F)]*aucMI
length(which(MI!=0))


MIwei <- cbind(colnames(data), MI)
# write.csv(MIwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/MIwei.csv", row.names = F)
MIfea <- colnames(data)[which(MI!=0)]
# write.csv(MIfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/MIfea.csv", row.names = F)



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
Relief <- imprelief[order(rankrelief, decreasing = F)]


## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(Relief == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(Relief == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucRelief <- picauc$auc

Relief <- imprelief[order(rankrelief, decreasing = F)]*aucRelief
length(which(Relief!=0))


Reliefwei <- cbind(colnames(data), Relief)
# write.csv(Reliefwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/Reliefwei.csv", row.names = F)
Relieffea <- colnames(data)[which(Relief!=0)]
# write.csv(Relieffea, file = "D:/E/博士/R_程序/EFS/Results/Feature/Relieffea.csv", row.names = F)



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

DT <-  impdt[order(rankdt, decreasing = F)]
# View(DT)

## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(DT == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(DT == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucDT <- picauc$auc

DT <- impdt[order(rankdt, decreasing = F)]*aucDT
length(which(DT!=0))


DTwei <- cbind(colnames(data), DT)
# write.csv(DTwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/DTwei.csv", row.names = F)
DTfea <- colnames(data)[which(DT!=0)]
# write.csv(DTfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/DTfea.csv", row.names = F)


# mRMR --------------------------------------------------------------------

library(mRMRe)

setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA') 

## data
data_all <- read.table("TCGA_pro_outcome_TN_log_train.txt", header = T, check.names = FALSE, sep="\t")
# View(data_all[,1:10])
dim(data_all)    # 158 824
Data <- data.frame(data_all)
# View(data[,1:10])


Data[,1] <- factor(Data[,1])
Data[,1] <- ordered(Data[,1], levels = c("0", "1"))


x <- Data[,-1]    # sample*gene
# View(x[,1:10])
y <-  as.numeric(as.factor(Data[,1]))-1 

seed <- 123    # svmfs 函数必须的
data.mrmre.train <- mRMR.data(data=Data[,-1], strata = Data[,1])
## classical mRMR feature selection
# res.fs.mrmr <- mRMR.classic(data=data.mrmre.train, target_indices=1, feature_count=25)
## 与上等价：solution_count不是传统mRMR的参数，所以如果设置为1，它仅执行一个经典的mRMR
res.fs.mrmr <- mRMR.ensemble(data.mrmre.train, target_indices=1, feature_count=200, solution_count = 1) # , method="exhaustive"
list(apply(solutions(res.fs.mrmr)[[1]], 2, function(x, y) 
{ return(y[x]) }, y = featureNames(data.mrmre.train)))
# myfs <- c(myfs, list(apply(solutions(res.fs.mrmr)[[1]], 2, function(x, y) 
# { return(y[x]) }, y = featureNames(data.mrmre.train))))
selected.features.mrmre <- mRMRe::solutions(res.fs.mrmr)


## 没选出的为0
imp = cbind(1:norm, rep(0, norm))
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



## train and test
xtrain <- data.frame(datatrain[,-1])
xtrain[, which(mRMR == 0)] <- 0
# View(xtrain[,1:10])
ytrain <- datatrain[,1]
xtest <- data.frame(datatest[,-1])
xtest[, which(mRMR == 0)] <- 0
# View(xtest[,1:10])
ytest <- datatest[,1]


## fit and predict
# model <- glm(ytrain~., data = xtrain, family = "binomial")
model <- glm(ytrain~., data = xtrain, family = "binomial", control = list(maxit = 100))
summary(model)
p_test <- predict(model, xtest, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, ytest)
names(A_test)<- c("p", "outcome")
picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
aucmRMR <- picauc$auc

mRMR <- impmRMR[order(rankmRMR, decreasing = F)]*aucmRMR
length(which(mRMR!=0))


mRMRwei <- cbind(colnames(data), mRMR)
# write.csv(mRMRwei, file = "D:/E/博士/R_程序/EFS/Results/Weight/mRMRwei.csv", row.names = F)
mRMRfea <- colnames(data)[which(mRMR!=0)]
# write.csv(mRMRfea, file = "D:/E/博士/R_程序/EFS/Results/Feature/mRMRfea.csv", row.names = F)


# 整合 ----------------------------------------------------------------------


table=rbind(Median, Variance, Chi2, Relief, P_cor, S_cor, MI, mRMR, Ridge, DT, RF_Gini, RF_Acc)
colnames(table) <- colnames(data)
View(table[,1:10])
number <- dim(table)[1]
table = table/number
View(table[,1:10])




# Plot Figure -------------------------------------------------------------


name <- "BRCA"


## 將所有降序排列，选取30个genes
efs_table <- table
# efs_table <- table[,1:30]
order = TRUE

if(order == TRUE){
  efs_table <- efs_table[, order(colSums(efs_table), decreasing = T)]  # 降序
}
efs_table <- efs_table[,1:24]



## 原始的，对所有的进行生序排列
# efs_table <- table[,1:30]
# order = TRUE
# 
# if(order == TRUE){
#   efs_table <- efs_table[, order(colSums(efs_table))]  # 升序
# }



paranr = length(efs_table[1,])  
if(paranr>100){
  b =  colSums(efs_table)
  #b= a[order(a)]
  
  # pdf(paste(name,'.pdf', sep=""),
  #     width= 12,
  #     height= 12)
  barplot(b,
          ylim=c(0,1),
          main= 'Ensemble Feature Selection',
          xlab = "Features",
          ylab = "Importance values",
          axisnames = FALSE
  )
  # dev.off()
}

# View(b)
# setwd('D:\\E\\博士\\R_程序\\EFS\\Results')
# write.csv(b, file = "rank.csv")


if(paranr<35){h=10}else {h=(paranr/5)}

names=colnames(efs_table)
setwd('D:\\E\\博士\\R_程序\\EFS\\Results')
# write.csv(names, file = "bio.csv", row.names = F)

# cols = c('goldenrod1','SlateBlue',
#          'royalblue','indianred3',
#          'darkolivegreen1','darkgreen',
#          'MediumSeaGreen','DeepSkyBlue',
#          'HotPink', 'BlueViolet',
#          'Sienna2', 'LightYellow3')

cols = c('#73D6C7','#FFFEA9',
         '#BFB9DD','#FF776C',
         '#70B2D7','#FFB14F',
         '#A6E152','#FFCAE6',
         '#D9D9D9','#C87CC1',
         '#C3EDC1','#FFED53')

# pdf(paste(name,'.pdf', sep=""),
# width= 12,
# height= h)
par(mar=c(5, 4, 4, 10), xpd=TRUE)
barplot= barplot(efs_table,
                 xlim=c(0,1),
                 # main= 'Ensemble Feature Selection',
                 horiz=T,
                 las=1,   # 2
                 # names.arg=abbreviate(names),
                 col=cols)
# legend("topright", legend=row.names(efs_table), col=cols, lty=1, lwd=3 )
legend("bottomright", inset=c(0.05,0.12), legend=row.names(efs_table), col=cols, lty=1, lwd=3 )
text(colSums(efs_table)+0.035,barplot,
     format(round(colSums(efs_table), 3),T))
segments(1, 0, 1, 1.25*paranr, lty = 3, col = "gray40")
# dev.off()






