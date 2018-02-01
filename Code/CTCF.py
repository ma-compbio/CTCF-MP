from optparse import OptionParser
import pandas as pd
import numpy as np
import sys
import warnings

chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
	"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
	"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

threshold = 0

def CTCF_ChIP(cell,ChIP_name,select = False):
	peak = pd.read_table("../Data/%s/%s.bed" %(cell,ChIP_name),sep = "\t",header = None)
	if np.asarray(peak).shape[1] == 10:
		peak.columns = ['chromosome','start','end','name','score','strand','signalValue','pValue','qValue','summit']

	elif np.asarray(peak).shape[1] == 9:
		print "broad peak"
		peak.columns = ['chromosome','start','end','name','score','strand','signalValue','pValue','qValue']
		peak["summit"] = (peak['end'] - peak['start']) * 0.5
	peak["covered"] = 0
	peak["point"] = peak["start"] + peak["summit"]
	peak["new_start"] = np.min((peak["point"] - threshold,peak["start"]),axis = 0)
	peak["new_end"] = np.max((peak["point"] + threshold,peak["end"]),axis = 0)

	try:
		CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell),sep = ",")
	except:
		CTCF = pd.read_table("../Data/%s/fimo.csv" %(cell),sep = ",")
		CTCF = CTCF.sort_values(by=['chromosome','start'])
		CTCF["active"] = 0

	if select == False:
		CTCF = CTCF[CTCF["active"] > 0]
		CTCF.index = xrange(len(CTCF))

	CTCF[ChIP_name] = 0
	total = len(CTCF)
	count = 0

	for chr in chrom_list:
		peak_value = []
		temp_CTCF = CTCF[CTCF["chromosome"] == chr]
		temp_peak = peak[peak["chromosome"] == chr]

		for i in xrange(len(temp_CTCF)):
			if (count % 100 == 0) and (i != 0):
				print "Dealing with %s: %d of %d\r" %(ChIP_name,count,total),
				sys.stdout.flush()
			index = temp_CTCF.index[i]
			start = CTCF["start"][index]
			end = CTCF["end"][index]
			summit = (start + end)/2
			a = temp_peak[(temp_peak["new_start"] <= summit) & (temp_peak["new_end"] >= summit)]

			if(len(a) != 0):
				peak["covered"][a.index] = 1
				peak_value.append(np.max(peak['signalValue'][a.index]))
			else:
				peak_value.append(0)
			count += 1
		CTCF[ChIP_name][temp_CTCF.index] = peak_value
	print "\nMotif covering ChIP %d of %d \n" %(np.sum(peak["covered"]),len(peak))

	if select:
		CTCF["active"] = ((CTCF["active"] + CTCF[ChIP_name]))
		print np.sum(CTCF["active"] > 0)
	CTCF.to_csv("../Temp/%s/CTCF.csv" %(cell),index=False)

def CTCF_CH(cell):
	loop = pd.read_table("../Data/%s/loop_from_Ch.csv" %(cell),sep = ",")

	chrom = np.concatenate((loop["chromosome1"],loop["chromosome2"]),axis = 0)
	start = np.concatenate((loop["start1"],loop["start2"]),axis = 0)
	end = np.concatenate((loop["end1"],loop["end2"]),axis = 0)
	start -= threshold
	end += threshold
	peak = {"chromosome":chrom,"start":start,"end":end}
	peak = pd.DataFrame(peak)
	peak["covered"] = 0

	CTCF = pd.read_table("../Data/%s/fimo.csv" %(cell),sep = ",")
	CTCF = CTCF.sort_values(by=['chromosome','start'])
	CTCF["active"] = 0
	#CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell),sep = ",")
	total = len(CTCF)
	count = 0

	for chr in chrom_list:
		active = []
		peak_value = []
		temp_CTCF = CTCF[CTCF["chromosome"] == chr]
		temp_peak = peak[peak["chromosome"] == chr]

		for i in xrange(len(temp_CTCF)):
			if (count % 100 == 0) and (i != 0):
				print "Generating Active Motif: %d of %d\r" %(count,total),
				sys.stdout.flush()
			index = temp_CTCF.index[i]
			start = CTCF["start"][index]
			end = CTCF["end"][index]
			summit = (start + end)/2
			a = temp_peak[(temp_peak["start"] <= summit) & (temp_peak["end"] >= summit)]

			if(len(a) != 0):
				active.append(1)
				peak["covered"][a.index] = 1
			else:
				active.append(0)
			count += 1
		CTCF["active"][temp_CTCF.index] = active
	print np.sum(CTCF["active"])
	CTCF.to_csv("../Temp/%s/CTCF.csv" %(cell),index=False)

def CTCF_Age(cell):
	CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell), sep = ",")
	#CTCF = CTCF[CTCF["active"] == 1]
	CTCF.index = xrange(len(CTCF))
	CTCF["age"] = 0
	Age = pd.read_table("../Data/%s/CTCF_age2" %(cell),sep = "\t")

	total = len(CTCF)

	age_dict={'hg19':1,'Human-Chimp':2,'Homininae':3,'Hominidae':4,
	'Catarrhini':5,'Simiiformes':6,'Haplorrhini':7,'Primate':8,
	'Strepsirrhini':9,'Human-Mouse':10,'Node_11':11,'Node_12':12,
	'Node_13':13,'Node_14':14,'Node_15':15,'Node_16':16,'Node_17':17,'Root':18,'NA':0,'nan':0}
	count = 0
	nomap = 0

	for chr in chrom_list:
		age_list = []
		temp_CTCF = CTCF[CTCF["chromosome"] == chr]
		temp_age = Age[Age["chromosome"] == chr]

		for i in xrange(len(temp_CTCF)):
			if (count % 100 == 0) and (i != 0):
				print "Mapping Age: %d of %d\r" %(count,total),
				sys.stdout.flush()
			index = temp_CTCF.index[i]
			start = CTCF["start"][index]
			end = CTCF["end"][index]
			summit = (start + end)/2
			a = temp_age[(temp_age["start"] <= summit) & (temp_age["end"] >= summit)]

			if(len(a) != 0):
				age_list.append(age_dict[str(a["age"][a.index[0]])])
			elif (len(a) > 1):
				print "error"
			else:
				nomap += 1
				age_list.append(0)
			count += 1
		CTCF["age"][temp_CTCF.index] = age_list

	print
	print "can't map: %d\n" %(nomap)


	CTCF["pos_index"] = 0 
	CTCF["neg_index"] = 0
	CTCF["pos_index"][CTCF["strand"] == "+"] = 1
	CTCF["neg_index"][CTCF["strand"] == "-"] = 1
	CTCF["pos_index"] = np.cumsum(CTCF["pos_index"])
	CTCF["neg_index"]  =np.cumsum(CTCF["neg_index"])

	CTCF.to_csv("../Temp/%s/CTCF.csv" %(cell),index = False)

def parse_args():
	parser = OptionParser(usage="CTCF Interaction Prediction", add_help_option=False)
	parser.add_option("-c","--cell",default = 'GM')
	
	(opts, args) = parser.parse_args()
	return opts

def run(cell):
	
	for CTCF in ['CTCF_peak']:
		CTCF_ChIP(cell,CTCF,True)
	
	for ChIP in ["DNase_Duke"]:
		CTCF_ChIP(cell,ChIP)
	
	CTCF_Age(cell)
	
def main():
	opts = parse_args()
	run(opts.cell)




if __name__ == '__main__':
	main()
