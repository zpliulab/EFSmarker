# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 22:53:43 2022

@author: LiLingyu
"""

import os
os.chdir('D:\E\博士\R_程序\EFS\Data\RTCGA')


##### 2022.2.16 ######
# https://machinelearningmastery.com/rfe-feature-selection-in-python/
# check scikit-learn version
import sklearn
print(sklearn.__version__)



#https://blog.csdn.net/weixin_43378396/article/details/90649064



import numpy as np
from pandas import DataFrame
#from sklearn.feature_selection import RFE
import pandas as pd


# load data
#dataframe = pd.read_table("TCGA_pro_outcome_TN_log_train.txt")  ## 79-0
dataframe = pd.read_table("TCGA_pro_outcome_TN_log_trainP.txt")
# sample name
featurenames = np.array(dataframe.columns)
# featire name
samplenames = np.array(dataframe.index)
# data form
esetarray = np.array(dataframe)
esetarray = esetarray.transpose()
# label
sampletype = [0]*79+[1]*79
# 
p = 2022


## Variance
from sklearn.feature_selection import VarianceThreshold
#sel = VarianceThreshold()
sel = VarianceThreshold(threshold=(2.5))
sel.fit_transform(dataframe)
varmatrix = sel.fit_transform(dataframe)
resultframe = DataFrame(varmatrix)
resultframe.to_csv("D:\E\博士\R_程序\EFS\DataPython/variance.txt", sep="\t")



## chi-square

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

model_sk = SelectKBest(score_func=chi2, k=50)
model_sk.fit(dataframe, sampletype)
print(model_sk.scores_)
print(model_sk.pvalues_)
ka2 = model_sk.scores_
ka2p = model_sk.pvalues_

resultframe = DataFrame(ka2)
resultframe.to_csv("D:\E\博士\R_程序\EFS\DataPython/chi2.txt", sep="\t")



## MI 


from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import mutual_info_classif


model_sk = SelectKBest(score_func=mutual_info_classif, k=50)
model_sk.fit(dataframe, sampletype)
print(model_sk.scores_)
mi = model_sk.scores_


resultframe = DataFrame(mi)
resultframe.to_csv("D:\E\博士\R_程序\EFS\DataPython/mi.txt", sep="\t")



## DT
## .ship[1]  -- 列数

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeClassifier

model_dtc = DecisionTreeClassifier()
kfold = KFold(n_splits=10, random_state=p)
scores = []

## 将dataframe转换成float的 esetarray
esetarray = np.array(dataframe)

for i in range(dataframe.shape[1]):
    score = cross_val_score(model_dtc, esetarray[:, i:i+1], sampletype, cv=kfold, scoring='mutual_info_score') # scoring='adjusted_mutual_info_score'
    scores.append((format(score.mean(), '.3f'), featurenames[i]))
print(scores)
dt = scores

resultframe = DataFrame(dt)
resultframe.to_csv("D:\E\博士\R_程序\EFS\DataPython/dt.txt", sep="\t")



## Relief


import numpy as np
from random import randrange
# from sklearn.datasets import make_classification
from sklearn.preprocessing import normalize

#三种范数的距离计算 
def distanceNorm(Norm, D_value):

	#距离范数
	if Norm == '1':
		counter = np.absolute(D_value)
		counter = np.sum(counter)
	elif Norm == '2':
		counter = np.power(D_value, 2)
		counter = np.sum(counter)
		counter = np.sqrt(counter)
	elif Norm == 'Infinity':
		counter = np.absolute(D_value)
		counter = np.max(counter)
	else:
		raise Exception('We will program this later......')
	return counter
 
#算法主体,features和lables格式为np.array,iter_ratio是迭代比例    
def fit(features, labels, iter_ratio):
	#初始化
	(n_samples, n_features) = np.shape(features)
	distance = np.zeros((n_samples, n_samples))
	weight = np.zeros(n_features)
 
	if iter_ratio >= 0.5:
		#计算距离
		for index_i in range(n_samples):
			for index_j in range(index_i+1, n_samples):
				D_value = features[index_i]-features[index_j]
				distance[index_i, index_j] = distanceNorm('2', D_value)
		distance += distance.T
	else:
		pass
  
	#开始迭代
	for iter_num in range(int(iter_ratio*n_samples)):
		#初始化
		nearHit = list()
		nearMiss = list()
		distance_sort = list()
 
		#随机抽取样本
		index_i = randrange(0, n_samples, 1)
		self_features = features[index_i]
 
		#搜索nearHit与nearMiss
		if iter_ratio >= 0.5:
			distance[index_i, index_i] = np.max(distance[index_i]) 
			for index in range(n_samples):
				distance_sort.append([distance[index_i, index], index, labels[index]])
		else:
			#分别计算距离
			distance = np.zeros(n_samples)
			for index_j in range(n_samples):
				D_value = features[index_i]-features[index_j]
				distance[index_j] = distanceNorm('2', D_value)
			distance[index_i] = np.max(distance)
			for index in range(n_samples):
				distance_sort.append([distance[index], index, labels[index]])
		distance_sort.sort(key = lambda x:x[0])
		for index in range(n_samples):
			if nearHit == [] and distance_sort[index][2] == labels[index_i]:
				nearHit = features[distance_sort[index][1]]
			elif nearMiss == [] and distance_sort[index][2] != labels[index_i]:
				nearMiss = features[distance_sort[index][1]]
			elif nearHit != [] and nearMiss != []:
				break
			else:
				continue
 
		#更新权重
		weight = weight-np.power(self_features-nearHit, 2)+np.power(self_features-nearMiss, 2)
	#print (weight/(iter_ratio*n_samples))
	return weight/(iter_ratio*n_samples)
 

#(features, labels) = make_classification(n_samples = 100)
## 将dataframe转换成float的 esetarray 
features = np.array(dataframe)    
labels = np.array( sampletype)
features = normalize(X = features, norm = 'l2', axis = 0)
#for x in range(1,10):
    #weight = fit(features,labels,0.8)
weight = fit(features, labels, 0.8)      
#if __name__ == '__main__':	test()

relief = weight

resultframe = DataFrame(relief)
resultframe.to_csv("D:\E\博士\R_程序\EFS\DataPython/relief.txt", sep="\t")



