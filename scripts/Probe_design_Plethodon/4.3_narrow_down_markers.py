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

# **** Requires numpy v1.8 *****

target_loci = 2430

filename = "4_spliced_transcripts.fa"
fh = open(filename,'r')
outfile = '4.3-reduced_final_seqs.fa'
fh_out = open(outfile,'a')

#intiate empty dictionary, we'll store the number of the marker as the key and all info from the list above as the value
exon_dict = {}
exon_list = []
filt_list = []
filt_dict = {}

#read in script 4 output of 2K + genes and write the seqs to a dict (exon_dict)
#I used this to capture the whole header, rather than SeqIO which just gets the seq name next to the >
for line in fh:
	line = line.strip()
	if line.startswith(">"):
		header = line
		exon_dict[header]=''
	else:
		seq = line
		exon_dict[header]=seq

#iterate through the dict and write the header and seq to a single line sep by "," and append to a list for the random sampling part		
for key,value in exon_dict.iteritems():
	if "LOC" not in key and "hypothetical" not in key and "putative" not in key:
		comb = key + "," + value
		exon_list.append(comb)


rand_list = random.sample(exon_list, target_loci)

#now write the parts of the final list to file
for item in rand_list:
	#convert to string to use .split
	line = str(item)
	#gene name
	gene = line.split(",")[0]
	#seq 
	seq = line.split(",")[1]
	#replace tabs with spaces
	gene = gene.replace("\t"," ")
	#write to output file
	fh_out.write(gene + "\n" + seq + '\n')

#count total number of loci
i = 0
for record in SeqIO.parse(outfile, "fasta"):
	i = i + len(record.seq)

print "total bp included in pared down target list = ", i

	
fh.close()
fh_out.close()



