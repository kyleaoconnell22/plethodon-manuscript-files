''''
Usage: python ID1Exon.py [input.fasta]

1) filters CSV file down to one exon per locus

Picks All exons from each gene, need to be filtered later if you only want one per gene

Filters for length so that we only keep the exons at least 200 bp
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

#open the gene list (Amby_gene_list.genePhred.txt)
filename = 'Amby_gene_list.genePhred.txt'
fh = open(filename, 'r')

#open the output file
file_out = '5_multi_exon_append.txt'
fh_out = open(file_out, 'a')

#filter the csv file (gene list)
for line in fh:
	if line.startswith('mRNA'):
		pass
	else:
		line = line.strip('')
		line = line.strip('\n')
		line = line.split('\t')
		gene = line[0]
		chrom = line[1]
		strand = line[2]
		exonS = line[8]
		exonS = exonS.split(',')
		#remove extra char at the end of both start and stop by excluding last char with -1
		exonS = exonS[0:-1]
		exonE = line[9]
		exonE = exonE.split(',')[0:-1]
		if len(exonS) == 1:
			if int(exonE[0])- int(exonS[0]) > 120:
				exonS = str(exonS[0])
				fh_out.write(gene + '_' + 'exon_1' + '\t' + chrom + '\t' + strand + '\t' + exonS + '\t' + str(exonE) +'\n')
		else:
			x = 0
			#iterate through the exon lists to make a new line for each exon, rather than listing them together or just picking one
			for i in range(0, len(exonS)):
				genesplit = gene.split('_')
				if len(genesplit) < 3:
					print genesplit
				else:
					
					gene1 = genesplit[0]
					gene2 = genesplit[1]
					gene3 = genesplit[2]
					genecomb = gene1 + '_' + gene2 + '_' + gene3
					gene_ex = genecomb + '_' + 'exon' + '_' + str(x)
					x = x + 1
					
					if int(exonE[i])-int(exonS[i]) > 120: #filter for length
						#write the output
						fh_out.write(gene_ex + '\t' + chrom + '\t' + strand + '\t' + exonS[i] + '\t' + exonE[i]+'\n')
			
		


#close files		
fh.close()
fh_out.close()


