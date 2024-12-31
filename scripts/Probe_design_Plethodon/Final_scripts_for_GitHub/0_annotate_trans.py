''''
Usage: python 0_annotate_trans.py [input.fasta]

assumes you have a text file with gene annotations and then appends these annotations to the transcript file

Requires single line transcript

Skip this if your transcript is already annotated
------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
Sept 2018
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

#assign base dir and change to it
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#open the annotation file
filename = "sals.1fpkm_blast_descs.txt" #this is my text file with a single line of annotation notes
'''
looks like this
TR8|c0_g1	macrophage-expressed gene 1 protein
TR63|c0_g1	tetratricopeptide repeat protein 26
TR74|c2_g1	reverse transcriptase-like protein
TR90|c0_g3	rna-directed dna polymerase from mobile el
'''
fh = open(filename, 'r')

#open the original transcriptome
transcriptome = "Sals.1fpkm.5ind_singleline.fasta"

#open the annotated outfile transcriptome
file_out = 'Sals_annotated.fasta'
'''
like this
>TR6|c2_g3
GCACCATAAACACTTACCCACCGAGCATGCACAGC
>TR6|c2_g4
GCACCATAAACACTTACCCACCGAGCATGCACAGC
'''
fh_out = open(file_out, 'a')

annodict = {}

#filter the csv file (gene list)
for line in fh: #TR8|c0_g1	macrophage-expressed gene 1 protein
	line = line.strip('')
	line = line.strip('\n')
	line = line.split('\t')
	trans = line[0]
	gene = line[1]
	gene = gene.replace(" ", "_")
	annodict[trans]=gene
fh.close()

i = 0
j = 0
#read transcriptome into SeqIO, and then write it to outfile with the annotation
for record in SeqIO.parse(transcriptome, "fasta"):
	for key, value in annodict.iteritems():
		if record.id == key:
			fh_out.write(">"+str(key)+"_"+str(value)+'\n'+str(record.seq)+"\n")
			i = i + 1

#print stats
print "{0} transcripts now have annotations".format(i)			

print "{0} transcripts have no annotations".format(j)	
fh.close()
fh_out.close()


