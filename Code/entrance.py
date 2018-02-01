import genInteraction_new
import genNegativeData
import genLabelData
import genVecs
import train
import warnings
import sys
import os
import CTCF
from optparse import OptionParser


def parse_args():
	parser = OptionParser(usage="CTCF Interaction Prediction", add_help_option=False)
	parser.add_option("-f", "--feature", default=100, help="Set the number of features of Word2Vec model")
	parser.add_option("-w","--word",default = 6)
	parser.add_option("-r","--range",default = 250)
	parser.add_option("-c","--cell",default = 'gm12878')
	parser.add_option("-t","--total",default = False)
	parser.add_option("-d","--direction",default = 'conv')


	(opts, args) = parser.parse_args()
	return opts

def makepath(cell,direction):
	if os.path.exists("../Temp/%s" %(cell)) == False:
		os.makedirs("../Temp/%s" %(cell))
	if os.path.exists("../Temp/%s/%s" %(cell,direction)) == False:
		os.makedirs("../Temp/%s/%s" %(cell,direction))

def run(word,feature,range,cell,total,direction):
	warnings.filterwarnings("ignore")
	makepath(cell,direction)
	if total!=False:
		print "Dealing with CTCF Motif Databese"
		CTCF.run(cell)

	if not os.path.isfile("../Temp/%s/CTCF.csv" %(cell)):
		print "Dealing with CTCF Motif Databese"
		CTCF.run(cell)

	if not os.path.isfile("../Temp/%s/CH.csv" %(cell)):
		print "Mapping Motifs to CHIA-PET"
		genInteraction_new.run(cell)

	if not os.path.isfile("../Temp/%s/%s/Negative.csv"%(cell,direction)):
		print "Generating Negative Data"		
		genNegativeData.run(cell,direction)
	
	if not os.path.isfile("../Temp/%s/%s/LabelData.csv"%(cell,direction)):
		genLabelData.run(cell,direction)

	if not os.path.isfile("../Temp/%s/Unsupervised"%(cell)):
		genVecs.Unsupervised(int(range),int(word),cell)
	if not os.path.isfile("../Temp/%s/%s/LabelSeq"%(cell,direction)):
		genVecs.gen_Seq(int(range),cell,direction)
	if not os.path.isfile("../Temp/%s/%s/datavecs.npy"%(cell,direction)):
		genVecs.run(word,feature,cell,direction)
	
	train.run(word,feature,cell,direction)
	
def main():
	opts = parse_args()
	run(opts.word,opts.feature,opts.range,opts.cell,opts.total,opts.direction)

if __name__ == '__main__':
	main()