
# 安装R包 --------------------------------------------------------------------

# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("curatedTCGAData")


# 加载R包 --------------------------------------------------------------------

library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# 设置路径 --------------------------------------------------------------------

setwd("D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA")


# 区分肿瘤样品和正常组织的表达量样品 -------------------------------------------------------

brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE) # 上四分位规范化RSEM TPM基因
brca <- TCGAutils::splitAssays(brca, c('01','11'))
xdata.raw <- t(cbind(assay(brca[[1]]), assay(brca[[2]])))

dim(xdata.raw)    # 1205 20501
# View(xdata.raw[,1:10])


# Get matches between survival and assay data
class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor
names(class.v) <- rownames(xdata.raw)

# View(rownames(xdata.raw))
# View( TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor )


# keep features with standard deviation > 0
# xdata_raw <- xdata.raw %>%
# { (apply(., 2, sd) != 0) } %>%
# { xdata.raw[, .] } %>%
#   scale
xdata_raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }

dim(xdata_raw)    # 1205 20222
# View(xdata_raw[,1:10])


# small.subset <- c('CD5', 'CSF2RB', 'HSF1', 'IRGC', 'LRRC37A6P', 'NEUROG2',
#                   'NLRC4', 'PDE11A', 'PIK3CB', 'QARS', 'RPGRIP1L', 'SDC1',
#                   'TMEM31', 'YME1L1', 'ZBTB11',
#                   sample(colnames(xdata.raw), 100))
# xdata <- xdata.raw[, small.subset[small.subset %in% colnames(xdata.raw)]]

small_subset <- colnames(xdata.raw)
dim(small_subset)    # 1205 20222
# View(small_subset[,1:10])

xdata <- xdata_raw[, small_subset[small_subset %in% colnames(xdata_raw)]]
dim(xdata)    # 1205 20222
# View(xdata[,1:10])
xdatat <- t(xdata)

ydata <- class.v
# View(ydata)
class(ydata)

# ydata[which(ydata == "Primary Solid Tumor")] <- c("Tumor")
# ydata[which(ydata == "Solid Tissue Normal")] <- c("Normal")

# save(xdata,ydata,file = 'TCGA_pro_outcome.rdata')


# 保存所有样本的分类数据 -------------------------------------------------------------

# xydata <- cbind(ydata, xdata)
# dim(xydata)    # 1205 20223
# View(xydata[,1:10])
# 
# xydata[which(xydata[,1] == "1"), 1] <- c("Tumor")
# xydata[which(xydata[,1] == "2"), 1] <- c("Normal")
# View(xydata[,1:10])

# write.table(xydata, "TCGA_pro_outcome.txt",quote=F,sep="\t")



# 112个Normal样本匹配和同一个体的Tumor样本进行匹配 -----------------------------------------

# which(ydata[,1] == "Primary Solid Tumor")


## 正常样本
sample_N <- rownames(xdata)[which(ydata == "Solid Tissue Normal")]
sample_N1 <- sample_N %>% 
  as_tibble() %>% 
  mutate(sample_N = substr(sample_N, 1, 12)) 
## 所有样本
sample <- rownames(xdata)
sample1 <- sample %>% 
  as_tibble() %>% 
  mutate(sample = substr(sample, 1, 12)) 

## 提取Tumor样本
lab <- which(as.matrix(sample1[,2]) %in% as.matrix(sample_N1[,2]))
xdata_TN <- xdata[lab,]
dim(xdata_TN)    # 224 20222
# View(xdata_TN[,1:10])

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")
# View(ydata_TN)

data_TN <- cbind(ydata_TN, xdata_TN)
dim(data_TN)    # 224 20223
# View(data_TN[,1:10])
# write.table(t(data_TN), "TCGA_pro_outcome_TN.txt",quote=F,sep="\t")




setwd("D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA")

# 加载R包 --------------------------------------------------------------------


# BiocManager::install("DESeq2")
library(DESeq2)
library(limma)
library(pasilla)



# 加载数据 --------------------------------------------------------------------

data <- read.table("TCGA_pro_outcome_TN.txt",header=T,sep='\t', check.names = F)
dim(data)   # 20223   224
# View(data[,1:10])
## 提取所有gene
# allgene <- rownames(data)[-1]
# View(allgene)
# write.csv(allgene, file = "allgene_list.csv", row.names = F)

xdatat <- t(as.matrix(apply(as.matrix(data[-1,]),1,function(x) as.numeric(x))))
colnames(xdatat) <- colnames(data)
# View(xdatat[,1:10])
xdatat[1,2]
xdata[2,1]
ydata <- data[1,]


class(xdatat)
class(xdata)

View(xdatat[,1:10])

# DESeq2包来对RNA-seq数据做差异分析 -------------------------------------------------

exprSet <- round(xdatat)
dim(exprSet)    # 20222  1205    20222   224
# View(exprSet[,1:10])

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")
# View(ydata_TN)

group_list <- as.factor(ydata_TN)
# group_list <- xydata[,1]
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)

dds2 <- DESeq(dds)  ##第二步，直接用DESeq函数即可
resultsNames(dds2)


# 提取你想要的差异分析结果，我们这里是Tumor组对Normal组进行比较
res <-  results(dds2, contrast=c("outcome","Tumor","Normal"))
summary(res) 
# plotMA(res)
# 
# # BiocManager::install("apeglm")
# library(apeglm)
# resLFC <- lfcShrink(dds2, coef="outcome_Tumor_vs_Normal", type="apeglm") #经过lfcShrink 收缩log2 fold change
# plotMA(resLFC) 

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_TN.csv")


# 确定阈值，筛选差异表达基因 -----------------------------------------------------------

## FDR较正
res1 <-  results(dds2, alpha = 0.01)  #默认FDR小于0.1，现取阈值padj小于0.05
# write.csv(res1,file= "DEG_res.csv")

diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > 3.321928) #
# log2(9)=3.169925,10-3.321928
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)

# write.csv(diff_gene_deseq2,file= "DEG_Tumor_vs_Normal_10.csv")





label <- c("BRCA1", "BRCA2", "PALB2", "PTEN", "TP53", 
           "ATM", "CDH1", "CHEK2", "NBN", "NF1", "STK11",
           "BARD1", "BRIP1", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM",
           "RAD51C", "RAD51D")

which(label %in% rownames(diff_gene_deseq2))    # 1  2 12 13 18




# 数据归一化 -------------------------------------------------------------------

# normalized_counts <- counts(dds2, normalized=TRUE)
# View(normalized_counts[,1:10])
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20222   224
View(dat111[,1:10])
# write.csv(dat111,file= "deseq_nor.csv")


## 标签是tumor和normal
data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
data0[3,3]


ydata_TN <- rbind(as.matrix(rep(c("0"), 112)), as.matrix(rep(c("1"), 112)))
colnames(ydata_TN) <- c("outcome")
data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
View(data0[,1:10])
# write.table(data0,"TCGA_pro_outcome_TN_log.txt",quote=F,sep="\t") 



# DESeq2软件独有的normlization方法 -----------------------------------------

# rld <- rlogTransformation(dds2)  ## 得到经过DESeq2软件normlization的表达矩阵！
# exprSet_new=assay(rld)
# par(cex = 0.7)
# n.sample=ncol(exprSet)
# if(n.sample>40) par(cex = 0.5)
# cols <- rainbow(n.sample*1.2)
# par(mfrow=c(2,2))
# boxplot(exprSet, col = cols,main="expression value",las=2)
# boxplot(exprSet_new, col = cols,main="expression value",las=2)
# hist(exprSet)
# hist(exprSet_new)



# 提取差异gene symbol ---------------------------------------------------------

Data <- read.table("TCGA_pro_outcome_TN_log.txt",header=T,sep='\t', check.names = F)
dim(Data)   # 20223   224
# View(Data[,1:10])
gene <- as.vector(rownames(Data)[-1])    # 20222
View(gene)    


DE_outcome <- read.csv("DEG_Tumor_vs_Normal_10.csv", header = T, sep=',')    # 489
Degene <- DE_outcome[,1] # 489
# write.csv(Degene, file = "DEgene_list.csv", row.names = F)


intersect(gene, Degene)


data <- Data[Degene,]
dim(data)    # 1633 1080
View(data[,1:10])


all_data <- rbind(Data[1,], data)
dim(all_data)    # 2568 1080
# View(all_data[,1:10])
# write.table(all_data, file = "TCGA_pro_outcome_TN_log_DE489.txt",quote = F, sep = "\t")


