import os
import sys

'''
rename dir of fasta files to locus names
secapr names them numbers and writes a key to reference_fasta_header_info.txt
This script just names them back
The reference fasta header needs to be in the list dir or give it the full path from s4
python 1_rename_fastas.py
Kyle O'Connell 
June 2020s
'''

fh_in = open('reference_fasta_header_info.txt','r')
header = []
for line in fh_in:
	line=line.strip()
	header.append(line)
	
for fasta in os.listdir('.'):
	if fasta.endswith('.fasta'):
		for locus in header:
			if fasta.split('.')[0] == locus.split('\t')[0]:
				newname = locus.split('\t')[1].split('|')[0].replace('/','-')+'.fasta'
				print fasta, newname
				os.rename(fasta,newname)