#!/usr/bin/python

import argparse
import os
import subprocess as sp
import shutil
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np


'''
This script will add missing taxa to alignments missing all individuals

uses seqio parse, so has extra part to rewrite in single line fasta

Kyle O'Connell
kyleaoconnell22@gmail.com
4/19/19

python add_missing_sequences.py -i . -o glocks_allsamples 

*New version removes sequences that have less than half of the bp because they often mess up the snp calling
*Also replaces - with N for SNP calling

'''
sample_set = []
missing_cutoff = 0.5

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: Specify the output dir to move the modified files")
    
    return parser.parse_args()

def open_input(in_dir):
	#get list of samples
	sample_list = []
	os.chdir(in_dir)
	for file in os.listdir('.'):
		for record in SeqIO.parse(file, "fasta"):
			sample_list.append(record.id)
	sample_set = set(sample_list)
	#sanity check
	print 'there are {0} phased samples in these alignments'.format(len(sample_set)/2)
	print '\n'

	#now do stuff
	for file in os.listdir('.'):
		#create and open outfile
		outtemp = file.split('.')[0] + '.temp.fasta'
		fh_temp = open(outtemp, 'a')
		#create list to figure out if seq is in this particular alignment
		name_list = []
		seqs = []
		i = 0
		#tdict={}
		exclude = []
		#print good sequences, basically just copying the other file
		for record in SeqIO.parse(file, "fasta"):
			#tdict[record.id]=str(record.seq)
			
			#but remove seqs that have less than half of the sequence, which will remove some crappy ones
			miss = float(len(str(record.seq.strip('-'))))
			full = float(len(str(record.seq)))
			if miss < missing_cutoff:
				exclude.append(record.id)
				print miss
			else:
				name_list.append(record.id)
				seqs.append(len(record.seq))
		#make list of single alleles		
		for sample in sample_set:
			base_0 = sample.split('_')[0]+'_'+sample.split('_')[1]+'_'+sample.split('_')[2]+'_'+sample.split('_')[3]+'_0'
			base_1 = sample.split('_')[0]+'_'+sample.split('_')[1]+'_'+sample.split('_')[2]+'_'+sample.split('_')[3]+'_1'
			if base_0 not in name_list and base_1 not in name_list:
				pass
			elif base_1 not in name_list and base_0 in name_list:
				exclude.append(base_0)
			elif base_0 not in name_list and base_1 in name_list:
				exclude.append(base_1)
		for record in SeqIO.parse(file, "fasta"):
			if not record.id in exclude:
				SeqIO.write(record, fh_temp, "fasta")

				
				
#SeqIO writes multiline fastas, so this module just rewrites the file as single line
#but it introduces an extra newline at the top, so the module below removes that
def single_line(in_dir):
	os.chdir(in_dir)
	for file in os.listdir('.'):
		if file.endswith('.temp.fasta'):
			fh = open(file, 'r')
			file_out = file.split('.')[0] + '_final.fasta'
			fh_out = open(file_out,'a')
			for line in fh:
				#write header line
				if line.startswith('>'):
					fh_out.write('\n' + line)
				else:
					#write the sequence, but also replace the '-' with 'N'
					line = line.strip()
					line = line.replace('-','N')
					fh_out.write(line)

def move_outfiles(out_dir, in_dir):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	os.chdir(in_dir)
	for file in os.listdir('.'):
		if file.endswith('_final.fasta'):
			shutil.move(file, out_dir)	
		elif file.endswith('.temp.fasta'):
			os.remove(file)								

#The single_line script adds an extra white line, so this module removes that and writes a 
#new file called 'single'		
def strip_white_space(out_dir):
	os.chdir(out_dir)
	for file in os.listdir('.'):
		if file.endswith('_final.fasta'):
			fh = open(file, 'rw')
			fileout = file.split('.')[0]
			fileout = fileout.strip('_final')
			fileout = fileout + '_single.fasta'
			fh_out = open(fileout, 'a')
			for line in fh:
				if line.startswith('\n'):
					pass
				else:
					fh_out.write(line)
			os.remove(file)						
					
def main():
	#define the arguments
	args = get_args()
	#change into the base dir with the secapr output alignments
	#execute searching and writing function
	open_input(args.in_dir)
	#convert to single line fasta
	single_line(args.in_dir)
	#move outfiles to out dir
	move_outfiles(args.out_dir, args.in_dir)
	strip_white_space(args.out_dir)
	
if __name__ == '__main__':
    main()