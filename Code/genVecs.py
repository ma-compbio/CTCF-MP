#encoding:utf-8
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
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
from gensim import corpora, models, similarities

def getkmervec(DNAdata1,DNAdata2,kmerdict):
	counter = 0
	DNAFeatureVecs = np.zeros((len(DNAdata1),2*len(kmerdict.keys())), dtype="float32")
	
	for DNA in DNAdata1:
		if counter % 1000 == 0:
			print "DNA %d of %d\r" % (counter, len(DNAdata1)),
			sys.stdout.flush()

		for word in DNA:
			DNAFeatureVecs[counter][kmerdict[word]] += 1
		counter += 1
	print
	
	counter = 0
	for DNA in DNAdata2:
		if counter % 1000 == 0:
			print "DNA %d of %d\r" % (counter, len(DNAdata2)),
			sys.stdout.flush()
		for word in DNA:
			DNAFeatureVecs[counter][kmerdict[word]+len(kmerdict.keys())] += 1
		counter += 1
	print

	return DNAFeatureVecs

def genkmerdict(k):
	vocab = ['A','C','T','G']
	dict1 = {}
	result = [""]
	for i in xrange(k):
		temp_result = []
		for word in result:
			if len(word) != i:
				print word 
				print len(word),i
				continue
			else:
				for letter in vocab:
					temp_word = word + letter
					temp_result.append(temp_word)

		result = temp_result

	count = 0
	for word in result:
		dict1[word] = count 
		count += 1

	return dict1

def Unsupervised(Range,word,cell):
	unlabelfile = open("../Temp/%s/Unsupervised" %(cell),"w")
	CTCF = pd.read_table("../Data/%s/fimo2.csv" %(cell),sep = ",")

	list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
			"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
			"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

	dict_pos={}
	total = len(CTCF)
	print total
	for i in xrange(total):
		if (i % 100 == 0) and (i != 0):
			print "number of unsupervised: %d of %d\r" %(i,total),
			sys.stdout.flush()

		start = CTCF["start"][i] - 1 - Range
		end = CTCF["end"][i] + Range
		chromosome = CTCF["chromosome"][i]
		if chromosome not in list:
			continue
		if dict_pos.has_key(chromosome):
			strs = dict_pos[chromosome]
		else:
			strs = processSeq.getString("../../Chromosome/" + chromosome + ".fasta")
			dict_pos[chromosome] = strs
		strand = CTCF["strand"][i]
		edstrs = strs[start:end]

		if strand == "-":
			edstrs = edstrs[::-1]
			edstrs = processSeq.get_reverse_str(edstrs)

		if "N" in edstrs:
			continue

		sentence = processSeq.DNA2Sentence(edstrs,word)
		unlabelfile.write(sentence+"\n")
	print

def gen_Seq(Range,cell,direction):
	print "Generating Seq..."
	table = pd.read_table("../Temp/%s/%s/LabelData.csv" %(cell,direction), sep=',')
	print len(table)
	table.drop_duplicates()
	print len(table)
	label_file = open("../Temp/%s/%s/LabelSeq" %(cell,direction), "w")

	label_file.write("label\t"
				  "seq1\t"
				  "strand1\t"
				  "seq2\t"
				  "strand2\t"
				  "chr\t"
				  "motif1\t"
				  "motif2\n")

	total = len(table)

	list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
			"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
			"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]
	number_positive = 0
	dict_pos={}

	cg1 = np.zeros((total),dtype="float32")
	cg2 = np.zeros((total),dtype="float32")
	cgm = np.zeros((total),dtype="float32")
	chr_length = np.zeros((total),dtype = "int64")
	for i in xrange(total):
		
		if (number_positive % 100 == 0) and (number_positive != 0):
			print "number of seq: %d of %d\r" %(number_positive,total),
			sys.stdout.flush()

		chromosome = table["chromosome"][i]
		if dict_pos.has_key(chromosome):
			strs = dict_pos[chromosome]
		else:
			strs = processSeq.getString("../../Chromosome/" + str(chromosome) + ".fasta")
			dict_pos[chromosome] = strs

		start1 = int(table["C1_start"][i] - 1 - Range)
		end1 = int(table["C1_end"][i] + Range)
		start2 = int(table["C2_start"][i] - 1 - Range)
		end2 = int(table["C2_end"][i] + Range)

		motif1 = table["C1_motif"][i]
		motif2 = table["C2_motif"][i]
		strand1 = table["C1_strand"][i]
		strand2 = table["C2_strand"][i]
		
		edstrs1 = strs[start1 : end1]
		edstrs2 = strs[start2 : end2]
		edstrs3 = strs[end1 : start2]

		
		chr_length[number_positive] = len(strs)

		if strand1 == "-":
			edstrs1 = edstrs1[::-1]
			edstrs1 = processSeq.get_reverse_str(edstrs1)
		if strand2 == "-":
			edstrs2 = edstrs2[::-1]
			edstrs2 = processSeq.get_reverse_str(edstrs2) 
		cgm[number_positive] = processSeq.countCG(edstrs3)
		cg1[number_positive] = processSeq.countCG(edstrs1)
		cg2[number_positive] = processSeq.countCG(edstrs2)
		
		if "N" in edstrs1 or "N" in edstrs2:
			table = table.drop(i)
			continue

		outstr = "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(table["label"][i],edstrs1,strand1,edstrs2,strand2,str(chromosome),motif1,motif2)
		label_file.write(outstr)
		number_positive += 1
	print
	table["CG1"] = cg1[:number_positive]
	table["CG2"] = cg2[:number_positive]
	table["CGM"] = cgm[:number_positive]
	print len(table)
	table.to_csv("../Temp/%s/%s/LabelData_select.csv" %(cell,direction),index=False)

#load word2vec model or train a new one
def getWord_model(word,num_features,min_count,cell):
	word_model = ""
	if not os.path.isfile("../Temp/model"):
		sentence = LineSentence("../Temp/%s/Unsupervised" %(cell),max_sentence_length = 15000)
		print "Start Training Word2Vec model..."
		# Set values for various parameters
		num_features = int(num_features)	  # Word vector dimensionality
		min_word_count = int(min_count)	  # Minimum word count
		num_workers = 20		 # Number of threads to run in parallel
		context = 20			# Context window size
		downsampling = 1e-3	 # Downsample setting for frequent words

		# Initialize and train the model
		print "Training Word2Vec model..."
		word_model = Word2Vec(sentence, workers=num_workers,\
						size=num_features, min_count=min_word_count, \
						window=context, sample=downsampling, seed=1,iter = 50)
		word_model.init_sims(replace=False)
		word_model.save("../Temp/model")
		#print word_model.most_similar("CATAGT")
	else:
		print "Loading Word2Vec model..."
		word_model = Word2Vec.load("../Temp/model")
		#word_model.init_sims(replace=True)


	return word_model

def getDNA_split(DNAdata,word):

	DNAlist1 = []
	DNAlist2 = []
	counter = 0
	for DNA in DNAdata["seq1"]:
		if counter % 100 == 0:
			print "DNA %d of %d\r" % (counter, 2*len(DNAdata)),
			sys.stdout.flush()

		DNA = str(DNA).upper()
		DNAlist1.append(processSeq.DNA2Sentence(DNA,word).split(" "))

		counter += 1

	for DNA in DNAdata["seq2"]:
		if counter % 100 == 0:
			print "DNA %d of %d\r" % (counter, 2*len(DNAdata)),
			sys.stdout.flush()

		DNA = str(DNA).upper()
		DNAlist2.append(processSeq.DNA2Sentence(DNA,word).split(" "))

		counter += 1
	print
	return DNAlist1,DNAlist2

def getAvgFeatureVecs(DNAdata1,DNAdata2,model,num_features, word,cell,direction):
	counter = 0
	DNAFeatureVecs = np.zeros((len(DNAdata1),2*num_features), dtype="float32")
	
	for DNA in DNAdata1:
		if counter % 1000 == 0:
			print "DNA %d of %d\r" % (counter, len(DNAdata1)),
			sys.stdout.flush()

		DNAFeatureVecs[counter][0:num_features] = np.mean(model[DNA],axis = 0)
		counter += 1
	print
	
	counter = 0
	for DNA in DNAdata2:
		if counter % 1000 == 0:
			print "DNA %d of %d\r" % (counter, len(DNAdata2)),
			sys.stdout.flush()
		DNAFeatureVecs[counter][num_features:2*num_features] = np.mean(model[DNA],axis = 0)
		counter += 1
	print

	return DNAFeatureVecs

def run(word, num_features,cell,direction):
	warnings.filterwarnings("ignore")

	global word_model,data

	word = int(word)
	num_features = int(num_features)
	word_model=""
	min_count=10

	word_model = getWord_model(word,num_features,min_count,cell)

	# Read data
	data = pd.read_table('../Temp/%s/%s/LabelSeq' %(cell,direction),sep = "\t")
	print "Generating Training and Testing Vector"
	
	
	datawords1,datawords2 = getDNA_split(data,word)
	dataDataVecs = getAvgFeatureVecs(datawords1,datawords2,word_model,num_features,word,cell,direction)

	print dataDataVecs.shape
	np.save("../Temp/%s/%s/datavecs.npy" %(cell,direction),dataDataVecs)
	
	'''
	dict1 = genkmerdict(6)
	datawords1,datawords2 = getDNA_split(data,6)
	dataDataVecs = getkmervec(datawords1,datawords2,dict1)
	np.save("../Temp/%s/%s/kmervecs.npy" %(cell,direction),dataDataVecs)
	'''

if __name__ == "__main__":
	run(6,100,'gm12878','conv')