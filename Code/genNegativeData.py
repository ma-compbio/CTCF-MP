import pandas as pd
import numpy as np
import warnings
import sys
import os

def gethead(head1,head2):
	head = []
	head.append('label')
	for name in head1:
		name = 'C1_'+name
		head.append(name)

	for name in head2:
		if name != "chromosome":
			name = 'C2_' + name
		head.append(name)
	return head

def getCTCF1(CTCF,index):
	temp = CTCF.iloc[index,0:]
	nonpredictors = ['name','chromosome','sth2','sth3','used']
	predictors_df = temp.drop(nonpredictors,axis = 1)
	head = list(predictors_df.columns.values)

	return np.asarray(predictors_df),head

def getCTCF2(CTCF,index):
	temp = CTCF.iloc[index,0:]
	nonpredictors = ['name','sth2','sth3','used']
	predictors_df = temp.drop(nonpredictors,axis = 1)
	head = list(predictors_df.columns.values)
	return np.asarray(predictors_df),head

def run(cell,direction):
	warnings.filterwarnings("ignore")
	CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell),sep = ",")
	Positive_data = pd.read_table("../Temp/%s/CH.csv" %(cell),sep = ",")

	if direction == 'conv':
		Positive_data = Positive_data[(Positive_data["C1_strand"] != "-") & (Positive_data["C2_strand"] != "+")]
	elif direction == 'tandem':
		Positive_data = Positive_data[Positive_data["C1_strand"] == Positive_data["C2_strand"]]
	elif direction == 'total':
		#Positive_data = Positive_data[(Positive_data["C1_strand"] != "-") | (Positive_data["C2_strand"] != "+")]
		Positive_data = Positive_data

	Positive_data = Positive_data[Positive_data["label"] == 1]
	Positive_data.index = xrange(len(Positive_data))

	distance = Positive_data["C2_start"] - Positive_data["C1_start"]

	distance = np.sort(distance,axis = None)

	#covering 95% percent of positive samples. Dump the 'extreme' situation
	a = int(len(Positive_data)*0.025)
	#a = 0
	min_range = distance[a]
	max_range = distance[len(distance)-a - 1]
	
	print min_range
	print max_range
	
	total = len(CTCF)
	CTCF.index = xrange(total)
	
	count = 0
	
	indexlist1=[]
	indexlist2=[]
	label = []

	if os.path.isfile("../Temp/%s/%s/neg_index.csv" %(cell,direction)):
		index = pd.read_table("../Temp/%s/%s/neg_index.csv" %(cell,direction),sep = ",")
		indexlist1 = index["index1"]
		indexlist2 = index["index2"]
		label = index["label"]
	else:

		for i in xrange(total):
			if (i % 100 == 0) and (i != 0):
					print "Processing: %d of %d\r" %(i,total),
					sys.stdout.flush()

			chromosome = CTCF["chromosome"][i]
			start = CTCF["start"][i]

			potential_pair = CTCF[(CTCF["chromosome"] == chromosome) & (CTCF["start"] <= start + max_range) & (CTCF["start"] >= start + min_range)]

			if len(potential_pair) != 0:
				for j in xrange(len(potential_pair)):
					indexlist1.append(i)
					indexlist2.append(potential_pair.index[j])
					if (CTCF["used"][i] == 1) & (CTCF["used"][potential_pair.index[j]] == 1):
						label.append(0)
					elif (CTCF["used"][i] != 1) & (CTCF["used"][potential_pair.index[j]] != 1):
						label.append(3)
					elif CTCF["used"][i] == 1:
						label.append(4)
					else:
						label.append(2)

					count +=1
		print
		print "Negative Data Source: " + str(count)

	indexlist1 = np.asarray(indexlist1)
	indexlist2 = np.asarray(indexlist2)
	label = np.asarray(label)
	label = label.reshape((len(label),1))
	a,head1 =  getCTCF1(CTCF,indexlist1)
	b,head2 = getCTCF2(CTCF,indexlist2)
	
	arrays = np.concatenate((label,a,b,np.abs(indexlist2-indexlist1).reshape((len(label),1))),axis = 1)
	header = gethead(head1,head2)
	header.append("num_between")
	table= pd.DataFrame(arrays,columns=header)
	table.to_csv("../Temp/%s/%s/Negative.csv" %(cell,direction), index = False)

	indexlist1 = indexlist1.reshape((len(indexlist1),1))
	indexlist2 = indexlist2.reshape((len(indexlist2),1))
	indexs = np.concatenate((label,indexlist1,indexlist2),axis = 1)
	indexs = pd.DataFrame(indexs,columns=['label','index1','index2'])
	indexs.to_csv("../Temp/%s/%s/neg_index.csv" %(cell,direction), index = False)
	
	CTCF.to_csv("../Temp/%s/CTCF.csv" %(cell),index = False)

if __name__ == "__main__":
	run()