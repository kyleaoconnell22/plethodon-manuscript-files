''''
Usage: python 2_script.py [1_single_exon.txt] [reference_genome.fa]

1) splices whole ref genome down to exons

------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
Sept 2018
------------------------
'''
#for this script to work as written, the genome has to be split by chromosome into sep. fasta files with a single seq in each
#Also, the dict method is really slow

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

#assign base dir and change to it, I recommend Trimmed
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#1_single_exon.txt
filename = sys.argv[2]
fh = open(filename, 'r')

#reference genome
filename2 = sys.argv[3]
fh2 = open(filename2, 'r')

#open the output file which is the spliced reference genomes
file_out = '2_spliced_reference.fa'
fh_out = open(file_out, 'a')
contig = []

#create empty dictionary where we will read in all sequences from the reference genome
Ref_dict = {}

#write the reference genome to a dict (so slow use a linux with lots of memory!)
for line in fh2:
	if line.startswith('>'):
		#split by the carrot to grab sample name
		line = line.split('>')
		name = line[1]
		contig.append(name)
		#add sample name as key to dict
		Ref_dict[name]=''
	else:
		#next line is sequence assign as dict value
		seq = line.upper()
		Ref_dict[name]+=seq

#print the number of contigs in the genome assembly
print "\n"
print "number of contigs in reference genome  = " , len(contig)
print "#####################################"
print "\n"

#make a variable to iterate at the end
#make two lists for counting at the end
i = 0
gene_list = []
chroms = []	

#for the single exon list (1_single_exon.txt)
for line in fh:
	#skip the header line
	if line.startswith('#') or line.startswith('Gene_Name'):
		pass
	#outer loop to assign variables
	else:
		line = line.strip('')
		line = line.split('\t')
		gene = line[0]
		chrom = line[1]
		strand = line[2]
		exonS = int(line[3])
		exonE = int(line[4])
		#lists for counting
		gene_list.append(gene)
		chroms.append(chrom)
		#iterate through the reference genome that was assigned to a dictionary
		for key in Ref_dict.keys():
			#if the header of the diction seq has the chrom (contig) name
			if key.startswith(chrom):
				#extract the exon from the full seq
				splice = Ref_dict[key]
				splice = splice[exonS-1:exonE-1]
				#i is the number of sequences we are splicing. It should match the number of genes at the end
				i = i + 1
				#write the pseudo reference genome
				fh_out.write('>' + gene + '|' + chrom + '\n' + splice + '\n')
				
chromset = set(chroms)
geneset = set(gene_list)		
print "number of unique contigs in gene list  = ", len(chromset)
print "number of (not unique) genes in gene list = ", len(gene_list) #note that this includes multiple exons per gene
print "number of unique genes in gene list = " , len(geneset)
print "number of sequences in new psuedoRef = " , i #This should be the total number of exons


fh.close()
fh2.close()
fh_out.close()