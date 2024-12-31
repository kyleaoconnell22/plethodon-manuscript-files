''''
Usage: python 1_IDexon_multi_genelist_same_names.py [input.csv]

filters CSV file down to one exon per locus

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
filename = sys.argv[2]
fh = open(filename, 'r')

#open the output file
file_out = '1_single_exon.txt'
fh_out = open(file_out, 'a')
#write the header
fh_out.write('#first exon in list in the case of mult. exons, if only one, then no choosing took place' + '\n')
fh_out.write('Gene_Name' + '\t' + 'Chromosome' + '\t' + 'Strand' + '\t' + 'Exon_Start'+'\t'+ 'Exon_END' + '\n')

#filter down the csv file
for line in fh:
	line = line.strip('')
	line = line.split('\t')
	gene = line[0]
	chrom = line[1]
	strand = line[2]
	exonS = line[8]
	exonS = exonS.split(',')
	exonS = exonS[0]
	exonE = line[9]
	exonE = exonE.split(',')
	exonE = exonE[0]
	fh_out.write(gene + '\t' + chrom + '\t' + strand + '\t' + exonS + '\t' + exonE + '\n')

#close files		
fh.close()
fh_out.close()


