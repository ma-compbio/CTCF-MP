#encoding:utf-8
import pandas as pd
import numpy as np
import os
import sys
import math
import random
import processSeq
import warnings

from sklearn import preprocessing
import sklearn.preprocessing
from sklearn.metrics import average_precision_score,precision_score,recall_score,f1_score
from sklearn.metrics import roc_auc_score,accuracy_score,matthews_corrcoef
from sklearn.cross_validation import cross_val_score,StratifiedKFold,train_test_split

import xgboost as xgb
from sklearn.externals import joblib
	
def analyzeResult(data,model,DataVecs):
	predict = model.predict(DataVecs)
	data['predict'] = predict
	data['predict_proba'] = model.predict_proba(DataVecs)[:,1]

	data['predict'][data['predict'] != 1] = 0
	data['label'][data['label'] != 1] = 0

	print ("Accuracy: %.4f %%" % (100. * sum(data["label"] == data["predict"]) / len(data["label"])))
	
	answer1 = data[data["label"] == 1]
	answer2 = data[data["label"] != 1]
	
	print ("Positive Accuracy: %.4f %%" % (100. * sum(answer1["label"] == answer1["predict"]) / len(answer1["label"])))
	print ("Negative Accuracy: %.4f %%" % (100. * sum(answer2["label"] == answer2["predict"]) / len(answer2["label"])))
	
	try:
		result_auc = model.predict_proba(DataVecs)[:,1]
		print ("Roc:%.4f\nAUPR:%.4f\n" % (roc_auc_score(data["label"],result_auc),
			average_precision_score(data["label"],result_auc)))
		print("Precision:%.4f\nRecall:%.4f\nF1score:%.4f\nMCC:%.4f\n" %(precision_score(data["label"],data["predict"]),
			recall_score(data["label"],data["predict"]),
			f1_score(data["label"],data["predict"]),
			matthews_corrcoef(data["label"],data["predict"])))
	except:
		print "ROC unavailable"

def score_func(estimator,X,Y):
	global accuracy,precision,recall,f1,mcc,auc,aupr,pos,neg,resultpredict,resultproba,resultlabel,count,mean_pos
	predict_proba = estimator.predict_proba(X)[:,1]
	True,False=1,0
	predict = estimator.predict(X)
	resultpredict = np.hstack((resultpredict,predict))
	resultproba = np.hstack((resultproba,predict_proba))
	resultlabel = np.hstack((resultlabel,Y))
	Y = (Y == 1)
	predict = (predict == 1)
	precision += precision_score(Y,predict)
	recall += recall_score(Y,predict)
	f1 += f1_score(Y,predict)
	accuracy += accuracy_score(Y,predict)
	mcc += matthews_corrcoef(Y,predict)
	auc += roc_auc_score(Y,predict_proba)
	aupr += average_precision_score(Y,predict_proba)
	pos += sum((predict == Y) & (Y == 1))/float(sum(Y==1))
	neg += sum((predict == Y) & ( Y== 0))/float(sum(Y==0))
	mean_pos += sum(predict_proba[Y == 1])
	count += 1
	print "finish %d of 10 accuracy: %.4f\tprecision: %.4f\trecall: %.4f\tf1: %.4f\tauc: %.4f\taupr: %.4f\r" %(count,accuracy/count,precision/count,recall/count,f1/count,auc/count,aupr/count),
	sys.stdout.flush()

	#joblib.dump(estimator, '../Models/boost_balance'+str(count)+'.model')

	return matthews_corrcoef(Y,predict)

def balance_data(data,dataDataVecs):
	np.random.seed(0)
	data.index = xrange(len(data))
	posdata = data[data["label"] == 1]
	negdata = data[data["label"] != 1]
	if len(posdata) > len(negdata):
		posdata = data[data["label"] != 1]
		negdata = data[data["label"] == 1]
	negdata["label"] = 0
	posdatavecs = dataDataVecs[posdata.index]
	negdatavecs = dataDataVecs[negdata.index]

	newnegindex = np.random.permutation(negdata.index)
	newnegindex = newnegindex[0:len(posdata)]
	negdata = negdata.reindex(newnegindex)
	negdatavecs = dataDataVecs[newnegindex]
	data = pd.concat([posdata,negdata])
	dataDataVecs = np.vstack((posdatavecs,negdatavecs))
	data.index = xrange(len(data))
	data["weight"] = 1
	return data,dataDataVecs

def setWeight(data):
	length = len(data)
	for i in xrange(5):
		a = data[data["label"] == i]
		if len(a) != 0:
			data["weight"][a.index] = length/float(len(a))		
	return data

def balance_data_dist(data,dataDataVecs):

	data["index_back"] = xrange(len(data))
	data["distance"] = data["C2_start"] - data["C1_start"]
	data_index = (data["distance"]/500).values
	data_index = np.floor(data_index)
	data["dist_window"] = data_index
	data["hash_value"] = 0
	count = 0

	positive = data[data["label"] == 1]
	positive = positive.sort_values(['distance'])
	negative = data[data["label"] !=1]
	negative = negative.sort_values(['distance'])
	positive.index = xrange(len(positive))
	negative.index = xrange(len(negative))
	hash_value = np.zeros((len(positive)),dtype="int64")
	for i in xrange(len(positive)):
		dist_window = positive["dist_window"][i]
		hash_value[i] = dist_window * 1000
		if i == 0:
			continue
		else:
			if positive["dist_window"][i] == positive["dist_window"][i-1]:
				count += 1
				hash_value[i] += count
			else:
				count = 0
	positive["hash_value"] = hash_value
	hash_value = np.zeros((len(negative)),dtype="int64")
	for i in xrange(len(negative)):
		hash_value[i] = negative["dist_window"][i] * 1000
		if i == 0:
			continue
		else:
			if negative["dist_window"][i] == negative["dist_window"][i-1]:
				count += 1
				hash_value[i] += count
			else:
				count = 0
	negative["hash_value"] = hash_value
	data = pd.concat([positive,negative])
	data = data.sort_values(["index_back"])
	data.index = xrange(len(data))
	positive = data[data["label"] == 1]
	negative = data[data["label"] != 1]
	positive.index = positive["hash_value"].values
	negative.index = negative["hash_value"].values
	positive = positive.drop_duplicates(["hash_value"])
	negative = negative.drop_duplicates(["hash_value"])

	positive_index = positive.index
	positive = positive.reindex(positive_index)
	negative = negative.reindex(positive_index)
	negative = negative.dropna()

	data = pd.concat([positive,negative])
	data.index = xrange(len(data))
	data["weight"] = 1
	print
	print len(positive)
	print len(negative)
	dataDataVecs = dataDataVecs[np.array(data["index_back"],dtype = 'int')]
	
	data.index = xrange(len(data))
	posdata = data[data["label"] == 1]
	negdata = data[data["label"] != 1]
	
	if len(posdata) > len(negdata):
		posdata = data[data["label"] != 1]
		negdata = data[data["label"] == 1]
	
	posdatavecs = dataDataVecs[posdata.index]
	negdatavecs = dataDataVecs[negdata.index]
	np.random.seed(0)
	newnegindex = np.random.permutation(negdata.index)
	newnegindex = newnegindex[0:len(posdata)]
	negdata = negdata.reindex(newnegindex)
	negdatavecs = dataDataVecs[newnegindex]
	data = pd.concat([posdata,negdata])
	dataDataVecs = np.vstack((posdatavecs,negdatavecs))
	data.index = xrange(len(data))
	data["weight"] = 1
	
	return data,dataDataVecs

def getLabelData(cell,direction):
	data = pd.read_table('../Temp/%s/%s/LabelData_select.csv' %(cell,direction),sep = ",")
	data["weight"] = 1

	data["distance"] = data["C2_start"] - data["C1_start"]
	length = float(len(data))

	for i in xrange(5):
		a = data[data["label"] == i]
		
		if len(a) != 0:
			data["weight"][a.index] = length/float(len(a))
		
		print "Label %d : %d\n" %(i,len(a))

	print len(data)
	return data

def getDataVecs(data):
	data["strand1"] = (data['C1_strand'] == "+")
	data["strand2"] = (data['C2_strand'] == "+")
	data["num1"] = np.abs(data["C2_pos_index"] - data["C1_pos_index"])
	data["num2"] = np.abs(data["C2_neg_index"] - data["C1_neg_index"])

	motif_feature_list = ['CTCF_peak','DNase_Duke','age','sth1']

	list1 = ['distance','CGM','CG1','CG2','num1','num2','strand1','strand2']
	for feature in motif_feature_list:
		list1.append("C1_"+feature)
		list1.append("C2_"+feature)
	
	dataDataVecs = data[list1]
	header = list(dataDataVecs.columns.values)
	print header
	return np.asarray(dataDataVecs),header

def run(word, num_features,cell,direction):

	warnings.filterwarnings("ignore")

	global accuracy,precision,recall,f1,mcc,auc,aupr,pos,neg,count,mean_pos
	accuracy,precision,recall,f1,mcc,auc,aupr,pos,neg,count,mean_pos = 0,0,0,0,0,0,0,0,0,0,0
	global resultpredict,resultproba,resultlabel
	resultpredict = np.zeros((1),dtype = "float32")
	resultproba = np.zeros((1),dtype = "float32")
	resultlabel = np.zeros((1),dtype = "float32")
	resultindex = np.zeros((1),dtype="int")
	global importance_sum
		
	data = getLabelData(cell,direction)
	
	dataDataVecs = np.load("../Temp/%s/%s/datavecs.npy" %(cell,direction))
	data = setWeight(data)
	print np.max(data['weight'])
	print np.min(data['weight'])
	extra,header = getDataVecs(data)
	print extra[4,:]


	dataDataVecs = np.concatenate((extra,dataDataVecs),axis = 1)

	scalar = sklearn.preprocessing.StandardScaler()
	dataDataVecs = scalar.fit_transform(dataDataVecs)
	print dataDataVecs.shape
	
	label = np.asarray(data["label"])
	weight = np.asarray(data["weight"])

	forest = xgb.XGBClassifier(max_depth=12,learning_rate=0.01, n_estimators = 500,nthread=30)

	fold_num = 10
	
	cv = StratifiedKFold(y = label, n_folds = fold_num, shuffle = True, random_state = 0)
	
	scores= cross_val_score(forest, dataDataVecs, label,\
							scoring = score_func, cv = cv,\
	 						n_jobs = 1,fit_params={'sample_weight': weight})
	
	print
	print "Accuracy:%.4f %%\nPositive Accuracy:%.4f %%\nNegative Accuracy:%.4f %%\n\nPrecision:%.4f\nRecall:%.4f\nF1:%.4f\nMCC:%.4f\n\nAUC:%.4f\nAUPR:%.4f\n " \
	%(accuracy*fold_num,pos*fold_num,neg*fold_num,precision/fold_num,recall/fold_num,f1/fold_num,mcc/fold_num,auc/fold_num,aupr/fold_num)