import sys
import os
import subprocess as sp
import shutil
import numpy
from Bio import SeqIO
'''
requires that magicblast is installed

'''
#give path to root
fasta_directory = sys.argv[1]
#move to the base directory
os.chdir(fasta_directory)

#give path to reference genome
ref_genome = sys.argv[2]
##########################################################
#create a blast db for the final fasta sheet
'''
proc_blastdb = sp.Popen(['makeblastdb', '-in', ref_genome, '-dbtype', 'nucl'])
proc_blastdb.wait()
'''
db = ref_genome



for filetype in os.listdir('.'):
	if filetype.endswith('_R1.fastq'):
		#print filetype
		query = filetype
		outfile = query + '_'+ ref_genome + '.blastout'
		print 'blasting ', filetype
		proc_blastn = sp.Popen(['magicblast', '-db', db, '-query', query, '-out', outfile, '-outfmt', "tabular"])
		proc_blastn.wait()
		fh = open(outfile, 'a')
		fh.write(query + ''+ 'blast results')
		fh.close()
s