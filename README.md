# CTCF-MP

## Required Package
CTCF-MP requires:
* Python (tested 2.7.13)
* numpy (tested 1.13.1)
* pandas (tested 0.20.3)
* gensim (tested 3.2.0)
* sklearn (tested 0.19.1)
* xgboost (tested 0.7)

## Required Data
Put hg19 sequences in a folder named 'Chromosome' and named each file as ‘chr1.fasta’,'chr2.fasta'...
The folder directory should look like this:

\Chromosome
\CTCF-MP
	\CTCF-MP\Code
	\CTCF-MP\Data
...


## Usage
The parameters are as followed:
* "-c","--cell",default = 'gm12878' :Controls the data it runs
* "-w", "--word",default = 6 :Controls k in k-mer we choose
* "-r", "--range",default = 250 :Controls the size of flanking region of the CTCF motif. (The final length of DNA sequence would be 2r+length of CTCF motif)
* "-d", "--direnction",default = 'conv' :Controls the subset of the dataset we run. 'conv' for 'convergent','tandem' for 'in tandem','imb' for 'imbalance'.

##Example

python entrance.py -c gm12878 -d conv -r 250