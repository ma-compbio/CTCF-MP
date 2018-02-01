#encoding:utf-8
import pandas as pd
import numpy as np
import warnings
import sys
import random
def balance_data(data,k):
	np.random.seed(0)
	data.index = xrange(len(data))
	posdata = data[data["label"] == 1]
	negdata = data[data["label"] != 1]
	if len(posdata) > len(negdata):
		posdata = data[data["label"] != 1]
		negdata = data[data["label"] == 1]

	newnegindex = np.random.permutation(negdata.index)
	newnegindex = newnegindex[0:len(posdata)*k]
	negdata = negdata.reindex(newnegindex)
	data = pd.concat([posdata,negdata])
	data.index = xrange(len(data))
	return data

def balance_data_dist(data):

	data["index_back"] = xrange(len(data))
	data["distance"] = data["C2_start"] - data["C1_start"]
	data_index = (data["distance"]/200).values
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
	print
	print len(positive)
	print len(negative)
	
	data.index = xrange(len(data))
	posdata = data[data["label"] == 1]
	negdata = data[data["label"] != 1]
	
	if len(posdata) > len(negdata):
		posdata = data[data["label"] != 1]
		negdata = data[data["label"] == 1]
	
	newnegindex = np.random.permutation(negdata.index)
	newnegindex = newnegindex[0:len(posdata)]
	negdata = negdata.reindex(newnegindex)
	data = pd.concat([posdata,negdata])
	data.index = xrange(len(data))
	data = data.drop(['distance','dist_window','hash_value'],axis = 1)
	return data


def run(cell,direction):
	Positive_data = pd.read_table("../Temp/%s/CH.csv" %(cell),sep = ",")

	Negative_data = pd.read_table("../Temp/%s/%s/Negative.csv" %(cell,direction),sep = ",")


	if (direction == 'conv') or (direction == 'imb'):
		Positive_data = Positive_data[(Positive_data["C1_strand"] != "-") & (Positive_data["C2_strand"] != "+")]
		Negative_data = Negative_data[(Negative_data["C1_strand"] != "-") & (Negative_data["C2_strand"] != "+")]

	elif direction == 'tandem':
		Positive_data = Positive_data[Positive_data["C1_strand"] == Positive_data["C2_strand"]]
		Negative_data = Negative_data[Negative_data["C1_strand"] == Negative_data["C2_strand"]]
	elif direction == 'total':
		#Positive_data = Positive_data[(Positive_data["C1_strand"] != "-") | (Positive_data["C2_strand"] != "+")]
		#Negative_data = Negative_data[(Negative_data["C1_strand"] != "-") | (Negative_data["C2_strand"] != "+")]
		Positive_data = Positive_data
		Negative_data = Negative_data

	table = pd.concat([Positive_data,Negative_data])
	table.index = xrange(len(table))
	table = table.drop_duplicates(["C1_start","C2_start","chromosome"],keep = False)
	Positive_data = Positive_data[Positive_data["label"] == 1]
	table = pd.concat([table,Positive_data])
	table.index = xrange(len(table))
	table = table.drop_duplicates(keep = "first")
	table = table[table["label"] != -1]
	if direction != "imb":
		table = balance_data_dist(table)
	else:
		table = table


	table.to_csv("../Temp/%s/%s/LabelData.csv" %(cell,direction),index=False)

	table = pd.read_table("../Temp/%s/%s/LabelData.csv" %(cell,direction),sep = ",")
	
	print len(table)

	print "Positive : %d Negative : %d\n" %(len(table[table["label"] == 1]),len(table[table["label"] != 1]))


if __name__ == "__main__":
	run()