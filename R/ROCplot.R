## clear
rm(list = ls())

## package
library(pROC)
library(ggplot2)
library(ROCR)  
## data
setwd('D:\\E\\博士\\R_程序\\EFS\\Data\\FeatureAUC')

roc1 <- read.csv(" feature_GSE21422 .csv", header=TRUE, sep = ',')
roc2 <- read.csv(" feature_GSE26910 .csv", header=TRUE, sep = ',')
roc3 <- read.csv(" feature_GSE45827 .csv", header=TRUE, sep = ',')
roc4 <- read.csv(" feature_GSE10797 .csv", header=TRUE, sep = ',')
roc5 <- read.csv(" feature_GSE10810 .csv", header=TRUE, sep = ',')
# roc6 <- read.csv(" feature_GSE65194 .csv", header=TRUE, sep = ',')
roc7 <- read.csv(" feature_GSE20437 .csv", header=TRUE, sep = ',')
roc8 <- read.csv(" feature_GSE61304 .csv", header=TRUE, sep = ',')
roc9 <- read.csv(" feature_GSE10780 .csv", header=TRUE, sep = ',')
roc10 <- read.csv(" feature_GSE42568 .csv", header=TRUE, sep = ',')
roc11 <- read.csv(" feature_GSE38959 .csv", header=TRUE, sep = ',')


myroc <- function(matrix){
  pred <- prediction(matrix[,1], matrix[,2])  
  perf <- performance(pred,"tpr","fpr") 
  x <- unlist(perf@x.values)  ##提取x值
  y <- unlist(perf@y.values)
  plotdata <- data.frame(x,y) 
  names(plotdata) <- c("x", "y")
  return(plotdata)
}

plotdata <- myroc(roc1)
ggplot(plotdata) + 
  geom_path(aes(x = x, y = y, colour = x), size=1) + 
  labs(x = "False positive rate", y = "Ture positive rate") +    # , title ="ROC曲线"
  scale_colour_gradient(name = 'False positive rate', low = 'blue', high = 'red') +
  # theme(plot.title = element_text(face = 'bold',size=15))
  theme_bw() + coord_equal()


Roc1 <- roc(roc1[,2],roc1[,1])
g <- ggroc(Roc1)
g
g + theme_minimal() + # ggtitle("输入标题") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="red", 
               linetype=6)###设置对角线，主要设置起点和终点的X,Y坐标
gl <- ggroc(Roc1, legacy.axes = TRUE)
gl

Roc2 <- roc(roc2[,2],roc2[,1])
Roc3 <- roc(roc3[,2],roc3[,1])
Roc4 <- roc(roc4[,2],roc4[,1])
Roc5 <- roc(roc5[,2],roc5[,1])
# Roc6 <- roc(roc6[,2],roc6[,1])
Roc7 <- roc(roc7[,2],roc7[,1])
Roc8 <- roc(roc8[,2],roc8[,1])
Roc9 <- roc(roc9[,2],roc9[,1])
Roc10 <- roc(roc10[,2],roc10[,1])
Roc11 <- roc(roc11[,2],roc11[,1])


g2 <- ggroc(list(GSE21422=Roc1, 
                 GSE26910=Roc2, 
                 GSE45827=Roc3, 
                 GSE10797=Roc4,
                 GSE10810=Roc5,
                 # GSE65194=Roc6, 
                 GSE20437=Roc7, 
                 GSE61304=Roc8,
                 GSE10780=Roc9,
                 GSE42568=Roc8,
                 GSE38959=Roc9),
            legacy.axes = TRUE)
g2


labels=c("GSE21422",
         "GSE26910",
         "GSE45827",
         "GSE10797",
         "GSE10810",
         # "GSE65194",
         "GSE20437",
         "GSE61304",
         "GSE10780",
         "GSE42568",
         "GSE38959")

g2 + annotate(geom = "segment", 
              x = 0, y = 0, xend =1, yend = 1, 
              colour = "gray", size = 0.5) +
  scale_fill_discrete(labels) +
  theme_gray() + coord_equal() +
  theme(legend.position = c(0.70,0.35), #legend.position = 'inside',
        legend.text = element_text(color = 'black',size = 10),
        axis.text = element_text(color = 'black',size = 15),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'))



# result0 <- read.csv("result_ROC.csv", header=TRUE, sep = ',')
result1 <- as.matrix(result0)
library(stargazer)
stargazer(result1)



