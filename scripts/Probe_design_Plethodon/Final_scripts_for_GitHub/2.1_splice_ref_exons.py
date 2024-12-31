'''
python 2.1_splice_ref_exons.py [Base_Dir] Notice that you need to hard code a bunch of file names

####### Edit this file name #######
means parameters that can be coded into the script

Edit everything within the parameters section below

Depending on which version of script 1 you used, it will splice the reference down to one exon per locus or all exons per locus. 

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

#assign base dir and change to it
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#1_single_exon.txt (note that this could be multi-exon-text)
filename = "1_single_exon.txt"
fh = open(filename, 'r')

#reference genome
#in this case the Ambystoma
reference_genome = 'AmexG_v3.0.0.fa'

#open the output file which is going to be the spliced reference genome
file_out = '2.1_spliced_reference_exons.fa'
fh_out = open(file_out, 'a')

#####

func = []

j = 0
i = 0
k = 0

exondict = {}

#parse the reference gene list to the functional genes we want
#AMEXTC_0340000293259_XYLT2_exon_1	AMEXG_0030000001	+	440201	742735	
for line in fh: #Full exon/gene list like the line above
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
		exondict[refgene]=line[0:5]
		k = k + 1
#print out the number of reference exons		
print "{0} possible exons from the reference gene list copied into the dictionary".format(k)

#open reference genome:
for record in SeqIO.parse(reference_genome, "fasta"):
	for key, value in exondict.iteritems():
		if record.id == value[1]:	
			exon = value[0]
			start = value[3]
			start = int(start)-1

			stop = value[4]
			stop = int(stop)-1
			
			#splice the reference genome
			splice = record.seq[start:stop]
			splice = splice.upper()
			seqlen = len(splice)
			if seqlen > 50:
				#write new spliced func to spliced ref genome outfile
				print "splicing and copying {0}".format(exon)
				fh_out.write('>' + str(exon) + '|' + str(record.id)+ ' ' + str(seqlen) + '\n' + str(splice) + '\n')

fh.close()
fh_out.close()
fh_mis.close()
print "###############################################"
print "###############################################"
print '\n'




