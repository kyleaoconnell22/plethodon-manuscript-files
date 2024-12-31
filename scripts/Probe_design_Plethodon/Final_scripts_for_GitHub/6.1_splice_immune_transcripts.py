''''
Usage: python ID1Exon.py [input.fasta]

assumes you have a text file with gene annotations and then appends these annotations to the transcript file

Requires single line transcript
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

#assign base dir and change to it, I recommend Trimmed
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#open the annotation file
filename = "sals.1fpkm_blast_descs.txt"
fh = open(filename, 'r')

#open the original transcriptome
transcriptome = "Sals.1fpkm.5ind_singleline.fasta"

#open the annotated outfile transcriptome
file_out = 'Sals_immunetranscripts.fasta'
fh_out = open(file_out, 'a')

immdict = {}
search_list = ["antigen", "mhc", "cytokine", "interleukin","chemokin","toll-like","antimicrobial peptide"]
#filter the csv file (gene list)
for line in fh: #TR8|c0_g1	macrophage-expressed gene 1 protein
	line = line.strip('')
	line = line.strip('\n')
	line = line.split('\t')
	trans = line[0]
	desc = line[1]
	desc = desc.replace(" ","_")
	for word in search_list:
		if word in desc:
			immdict[trans]=desc
i = 0
for record in SeqIO.parse(transcriptome, "fasta"):
	for key, value in immdict.iteritems():
		if record.id == key:
			fh_out.write(">"+str(key)+"_"+str(value)+'\n'+str(record.seq)+"\n")
			i = i + 1
print i
