'''
python 4_filter_blastout_splice_CDS [Base_Dir] Notice that you need to hard code a bunch of file names

####### Edit this file name #######
means parameters that can be coded into the script

Edit everything within the parameters section below
------------------------
written for Python 2.7
Blake O'connell and Kyle O'Connell
Sept. 2018
------------------------
'''
#________________________________________#
import sys
import os
import subprocess as sp
import random
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import numpy
import numpy as np
#________________________________________#
print '\n'
print "###############################################"
print "###############################################"

################!!!!!!!!!!!!!!!!!!Troubleshooting

#Parameters

fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#this is the functional gene list
filename = "vision_gene_list.txt"
fh = open(filename, 'r')

# open the second file which is the gene list, this can be all coding region or split into indiv exons using script 1_IDexon_multi_genelist
filehandle2 = '5_multi_exon_append.txt'
fh2 = open(filehandle2, 'r')

#assign reference genome for part below
reference_genome = "/Volumes/ANALYSIS/Plethodon_probes/Exon_ID_Pipeline/AmexG_v3.0.0.fa"

#assign output file
fileout = "functional_referencegenome.fa"
fh_out = open(fileout, 'a')

mismatch = "mismatch.out"
fh_mis = open(mismatch, 'a')
#####

func = []

j = 0

#write the functional gene list to a list
for line in fh: #for line in functional gene list
	line = line.strip()
	line = line.upper()
	func.append(line)
	j = j + 1

print "{0} functional genes are in your candidate list".format(j)

i = 0
k = 0
genedict = {}
funcdict = {}
funclist = []

keylist = []

#parse the reference gene list to the functional genes we want
#AMEXTC_0340000293259_XYLT2_exon_1	AMEXG_0030000001	+	440201	742735	
for line in fh2: #Full exon/gene list like the line above
	line = line.strip()
	#skip mRNA lines because they break the script below
	if line.startswith("mRNA"):
		pass
	#only use lines that start with the Ref gene, otherwise the script also breaks
	if line.startswith("AMEXTC"):
		line = line.split('\t')
		#AMEXTC_0340000023485_HOXB1_exon_1
		refgene = line[0]

		#create dict with key = target gene name
		#and value as an array of the full ref gene, contig loc, strand, start-stop of the coding seq
		genedict[refgene]=line[0:5]
		k = k + 1
		funclist.append(line)
print "{0} possible genes from the reference gene list copied into the dictionary".format(k)
#now see if the gene ID is in the functional gene list of target genes
#if so, write a new dictionary with just the coords of the functional genes 
countlist = []
for key, value in genedict.iteritems():
	key1 = key.split('_')[2]
	if key1 not in keylist:
		keylist.append(key1)
	if key1 in func:
		funcdict[key]=value
		i = i +1
		if key1 not in countlist:
			countlist.append(key1)
print "{0} reference genes matched exons in the functional gene list".format(i)
print "{0} reference genes matched genes in the functional gene list".format(len(countlist))


a = 0
for value in func:
	if value not in keylist:
		fh_mis.write(str(value) + '\n')
		a = a + 1

print "{0} reference genes mismatched the functional gene list and are printed to mismatch.out".format(a)

#open reference genome:
for record in SeqIO.parse(reference_genome, "fasta"):
	for key, value in funcdict.iteritems():
		if record.id == value[1]:			
			gene = value[0]
			start = value[3]
			start = int(start)-1
			
			stop = value[4]
			stop = int(stop)-1
			
			splice = record.seq[start:stop]
			seqlen = len(splice)
			otherlen = stop -start
			#write new spliced func ref genome
			#fh_out.write('>' + str(gene) + '|' + str(record.id)+ ' ' + str(seqlen) + '\n' + str(splice) + '\n')
fh.close()
fh_out.close()
fh_mis.close()
print "###############################################"
print "###############################################"
print '\n'




