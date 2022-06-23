
rm(list = ls())

# 数据输入 --------------------------------------------------------------------

setwd("D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA")

Data = read.table("TCGA_pro_outcome_TN_log.txt", header = T, check.names = FALSE)
dim(Data)    # 16145   224
# View(Data[,1:10])
gene <- as.vector(rownames(Data)[-c(1)])


DEData = read.table("TCGA_pro_outcome_TN_log_DE489.txt", header = T, check.names = FALSE)
dim(DEData)    # 391 1080
# View(DEData[,1:10])
Degene <-  as.vector(rownames(DEData)[-c(1)])

intersect(gene, Degene)  # 489
sub <- setdiff(gene, Degene)    # 从整体数据中，删除DEgene及其表达值
length(sub)

# biomarker ---------------------------------------------------------------

setwd('D:\\E\\博士\\R_程序\\EFS\\Data')

## mark
malacards <- as.matrix(read.csv("malacards82.csv", header = T, sep=','))

## scmark
# GEDFN <- as.matrix(read.csv("GEDFN169.csv", header = T, sep=','))

# 从 https://www.sciencedirect.com/topics/medicine-and-dentistry/tr --------
mama_70 <- as.matrix(read.csv("mamaprint70.csv", header = T, sep=',') )

KEGG_147 <- as.matrix(read.csv("KEGG147.csv", header = T, sep=',') )

tf <- as.matrix(read.csv("tfgene119.csv", header = T, sep = ','))

# 交集 ----------------------------------------------------------------------

intersect(malacards, tf)    # 2
intersect(malacards, mama_70)      # 0
intersect(malacards, KEGG_147)      # 10
intersect(mama_70, KEGG_147)    # 1
intersect(mama_70, tf)    # 0
intersect(KEGG_147, tf)    # 8


intersect(sub, malacards)      # 38
# intersect(sub, GEDFN)    # 165
intersect(sub, KEGG_147)    # 146
intersect(sub, mama_70)      # 52
intersect(sub, tf)      # 117


intersect(Degene, malacards)    # 0
# intersect(Degene, GEDFN)    # 3
intersect(Degene, mama_70)    # 4
intersect(Degene, KEGG_147)    # 1
intersect(Degene, tf)      # 0

# 102+9+47个gene整合 ------------------------------------------------------------

add_gene <- as.matrix( union( union(union(intersect(sub, malacards), 
                                          intersect(sub, mama_70)), 
                                    intersect(sub, KEGG_147)), 
                              intersect(sub, tf) ) 
                       )
length(add_gene)  # 334


rownames(add_gene) <- add_gene[,1] 

# View(Data[rownames(add_gene),1:5])    # 在wps中验证了

add_data <- Data[rownames(add_gene),]
dim(add_data)    # 412  224

all_data <- rbind(DEData, add_data)
dim(all_data)    # 824 224
View(all_data[,1:10])

geneall <- rownames(all_data)[-1]
View(geneall)

setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\RTCGA')
# write.csv(geneall, file = 'UNgene_list.csv', row.names = F)
# write.table(all_data, file = "TCGA_pro_outcome_TN_log_UN824.txt",quote = F, sep = "\t")


