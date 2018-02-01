#encoding:utf-8
import pandas as pd
import numpy as np 
import warnings
import sys
import os

def getcases(loop,i,CTCF,origin_CTCF):
	threshold = 0
	cases1 = ""
	i = loop.index[i]
	while (len(cases1) == 0) and (threshold <= 10):
		chromosome = loop["chromosome1"][i]
		start = loop["start1"][i]- threshold
		end = loop["end1"][i]+ threshold

		cases1 = CTCF[(CTCF["chromosome"] == chromosome) & (CTCF["start"] >= start - 18) & (CTCF["end"] <= end + 18)]
	
		if len(cases1) != 0:
			origin_CTCF["used"][cases1.index] = 1
		threshold +=50
	
	threshold = 0
	cases2 = ""
	while (len(cases2) == 0) and (threshold <= 10):

		chromosome = loop["chromosome2"][i]
		start = loop["start2"][i]- threshold
		end = loop["end2"][i]+ threshold

		cases2 = CTCF[(CTCF["chromosome"] == chromosome) & (CTCF["start"] >= start - 18) & (CTCF["end"] <= end + 18)]
		
		if len(cases2) != 0:
			origin_CTCF["used"][cases2.index] = 1
		threshold +=50
		

	return cases1,cases2

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

	print head
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

def run(cell):
	warnings.filterwarnings("ignore")
	loop = pd.read_table("../Data/%s/loop_from_Ch.csv" %(cell),sep = ",")
	CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell),sep = ",")

	total = len(loop)
	loop = loop[loop["chromosome1"] == loop["chromosome2"]]
	loop.index = xrange(len(loop))
	inter = total - len(loop)

	pos,unlab,upside,downside,unmap= 0,0,0,0,0
	indexlist1=[]
	indexlist2=[]
	label = []


	if os.path.isfile("../Temp/%s/pos_index.csv" %(cell)):
		index = pd.read_table("../Temp/%s/pos_index.csv" %(cell),sep = ",")
		indexlist1 = index["index1"]
		indexlist2 = index["index2"]
		label = index["label"]

	else:
		CTCF["used"] = 0
		chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
			"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
			"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

		count = 0
		for chr in chrom_list:
			temp_CTCF = CTCF[CTCF["chromosome"] == chr]
			temp_loop = loop[loop["chromosome1"] == chr]

			for i in xrange(len(temp_loop)):
				if (count % 100 == 0) and (i != 0):
					print "Mapping: %d of %d\r" %(count,len(loop)),
					sys.stdout.flush()
				count += 1

				cases1,cases2 = getcases(temp_loop,i,temp_CTCF,CTCF)
				length1,length2 = len(cases1),len(cases2)

				if length1 == 0 and length2 == 0:
					unmap += 1
					continue
				elif length1 == 0:
					upside +=1
					continue
				elif length2 == 0:
					downside +=1
					continue
				elif length1 == 1 and length2 == 1:
					label.append(1)
					indexlist1.append(cases1.index[0])
					indexlist2.append(cases2.index[0])
					pos += 1
				else:
					for j in xrange(length1):
						for k in xrange(length2):
							label.append(-1)
							indexlist1.append(cases1.index[j])
							indexlist2.append(cases2.index[k])
							unlab += 1

		print
		print "unmappable: %d  no upside: %d no downside : %d positive: %d potential: %d inter: %d\n" %(unmap,upside,downside,pos,unlab,inter)
	
	indexlist1 = np.asarray(indexlist1)
	indexlist2 = np.asarray(indexlist2)
	
	label = np.asarray(label)
	label = label.reshape((len(label),1))
	
	a,head1 =  getCTCF1(CTCF,indexlist1)
	b,head2 = getCTCF2(CTCF,indexlist2)

	indexlist1 = indexlist1.reshape((len(indexlist1),1))
	indexlist2 = indexlist2.reshape((len(indexlist2),1))
	indexs = np.concatenate((label,indexlist1,indexlist2),axis = 1)
	indexs = pd.DataFrame(indexs,columns=['label','index1','index2'])
	indexs.to_csv("../Temp/%s/pos_index.csv" %(cell), index = False)
	arrays = np.concatenate((label,a,b,np.abs(indexlist2-indexlist1)),axis = 1)
	header = gethead(head1,head2)
	header.append('num_between')
	table= pd.DataFrame(arrays,columns=header)
	table.to_csv("../Temp/%s/CH.csv" %(cell),index = False)

	
	loop.to_csv("../Data/%s/loop_from_Ch.csv" %(cell),index = False)
	CTCF.to_csv("../Temp/%s/CTCF.csv" %(cell),index = False)
if __name__ == "__main__":
	run('gm12878')