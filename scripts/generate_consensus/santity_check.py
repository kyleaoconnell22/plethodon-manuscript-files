import sys
import Bio
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
import os
import subprocess as sp
import shutil
import argparse
'''
Stand alone script for grabbing the sequence from a single locus for each species
requires that you run create_sp_consensus.py first
This is appended to the end of that script now, so don't need to run it sep
Kyle O'Connell
June 2020
'''
#module for testing the similarity between species sets
target = '348'
fh_out=open('sanity_check.fasta','a')
for filetype in os.listdir('.'):
	if filetype.endswith('consensus.fasta'):
		with open(filetype, "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				if record.id == target:
					fh_out.write('>'+filetype.split('_')[0]+'_')
					SeqIO.write(record,fh_out,'fasta')
fh_sanity = open('sanity_check.fasta','r')
fh2 = open(target+'_consensus.fasta','a')
for line in fh_sanity:
	if line.startswith('>'):
		line=line.strip()
		fh2.write('>'+line.split('_')[0].split('>')[1]+'_'+target+'\n')
	else:
		fh2.write(line)
os.remove('sanity_check.fasta')