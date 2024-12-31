#!/usr/bin/python

import argparse
import os
import subprocess as sp
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np


'''
This script will trim individuals missing more than the proportion defined as the missing_cutoff val
Also replaces - with N
Also filters alignments to designated samples so you can filter to one species for snp calling or, you can designate all samples
uses seqio parse, so has extra part to rewrite in single line fasta

keep file should be just a list
sample1
sample2
sample3

python trim_replacedash.py -i . -o glocks_allsamples -k keep_file.txt

Kyle O'Connell
kyleaoconnell22@gmail.com
4/19/19

'''
keep_list = []
sample_set = []
missing_cutoff = 0.50

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: Specify the output dir to move the modified files")
    parser.add_argument("-k", "--keep_file", required=True, help="REQUIRED: Specify the list of samples to keep")

    return parser.parse_args()

def keep_samples(keepfile):
	fh = open(keepfile, 'r')
	for line in fh:
		line = line.strip()
		keep_list.append(line)

def open_input(in_dir):
	#get list of samples
	sample_list = []
	os.chdir(in_dir)
	for file in os.listdir('.'):
		for record in SeqIO.parse(file, "fasta"):
			if not record.id.startswith('_R_'):
				sample_list.append(record.id)
	sample_set = set(sample_list)
	
	#sanity check
	print 'there are {0} samples in these alignments'.format(len(sample_set))
	print 'there are {0} samples in the keep file'.format(len(keep_list))
	print '\n'
	
	sample_out = 'sample_list.txt'
	fh_samp = open(sample_out,'w')
	for item in sample_set:
		 fh_samp.write(item+'\n')
	#now do stuff
	for file in os.listdir('.'):
		#create and open outfile
		outtemp = file.split('.')[0] + '.temp' + '.fasta'
		fh_temp = open(outtemp, 'a')
		#create list to figure out if seq is in this particular alignment
		name_list = []
		seqs = []
		i = 0
		
		#print 'printing good sequence'
		#print good sequences, basically just copying the other file
		for record in SeqIO.parse(file, "fasta"):
			#but remove seqs that have less than half of the sequence, which will remove some crappy ones
			if len(record.seq) > 50 and record.id in keep_list:
				miss = record.seq.count('-')
				miss = miss + 1
				perc = float(miss)/(len(record.seq))
			#round to two decimal places
				perc = round(perc,2)
				if perc > missing_cutoff:
					pass
				else:
					name_list.append(record.id)
					SeqIO.write(record, fh_temp, "fasta")
					seqs.append(len(record.seq))
			else:
				pass
		#set alignment length for this file	
		for item in seqs:
			i = item


#SeqIO writes multiline fastas, so this module just rewrites the file as single line
#but it introduces an extra newline at the top, so the module below removes that
def single_line(in_dir):
	print 'rewriting single lines'
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
	print 'moving outfiles'
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
	print 'removing white space'
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
	for file in os.listdir('.'):
		print 'renaming final files'
		if file.endswith('_single.fasta'):
			newname = file.split('_single.fasta')[0]+'.fasta'
			os.rename(file,newname)

def remove_empty(out_dir):
	os.chdir(out_dir)
	fh_empty=open('empty_files.txt','a')
	for fasta in os.listdir('.'):
		if fasta.endswith('.fasta'):
			print len(open(fasta).readlines())
			if len(open(fasta).readlines()) < 2:
				os.remove(fasta)
				fh_empty.write(fasta+'\n')
				print 'removing empty file ', fasta
def main():
	#define the arguments
	args = get_args()
	keep_samples(args.keep_file)
	#change into the base dir with the gblocks or trimal output alignments
	#execute searching and writing function
	open_input(args.in_dir)
	#convert to single line fasta
	single_line(args.in_dir)
	#move outfiles to out dir
	move_outfiles(args.out_dir, args.in_dir)
	strip_white_space(args.out_dir)
	remove_empty(args.out_dir)
if __name__ == '__main__':
    main()