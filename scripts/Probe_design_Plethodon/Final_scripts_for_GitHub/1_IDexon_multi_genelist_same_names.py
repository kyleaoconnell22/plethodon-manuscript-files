''''
Usage: python 1_IDexon_multi_genelist_same_names.py [input.csv]

filters CSV file down to one exon per locus

Picks All exons from each gene, need to be filtered later if you only want one per gene

The next script called 1_IDexon_single_genelist will just pick one exon per locus

Filters for length so that we only keep the exons at least 200 bp

Only parsing referene gene Phred file, not the actual reference genome
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
import shutil
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import numpy
import numpy as np
#________________________________________#

#assign base dir and change to it
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)
'''
this part assumes you are using the Ambystoma reference genome with the annotations in the genePhred format. If you use a GFF, all these scripts 
will need modification. Obviously, a different reference genome will require different names etc in some scripts
'''
#open the reference genome gene list (Amby_gene_list.genePhred.txt)
filename = sys.argv[2]
fh = open(filename, 'r')

#open the output file
file_out = '1_single_exon.txt'
fh_out = open(file_out, 'a')
#write the header
fh_out.write('Gene_Name' + '\t' + 'Chromosome' + '\t' + 'Strand' + '\t' + 'Exon_Start'+'\t'+ 'Exon_END' + '\n')

#filter the csv file (gene list)
for line in fh:
	line = line.strip('')
	line = line.strip('\n')
	line = line.split('\t')
	gene = line[0]
	chrom = line[1]
	strand = line[2]
	exonS = line[8]
	exonS = exonS.split(',')
	#remove extra char at the end
	exonS = exonS[0:-1]
	exonE = line[9]
	exonE = exonE.split(',')
	#iterate through the exon lists to make a new line for each exon, rather than listing them together or just picking one
	for i in range( 1, len(exonS)):
		if int(exonE[i])-int(exonS[i]) > 199:
			#write the output
			fh_out.write(gene + '\t' + chrom + '\t' + strand + '\t' + exonS[i] + '\t' + exonE[i]+'\n')
		
	

#close files		
fh.close()
fh_out.close()


