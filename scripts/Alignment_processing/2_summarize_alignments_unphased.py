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
Lots of this script was adapted from SuperCrunch by portik@github

Notice that
j = number of '-' char
i = total len of aln

Kyle O'Connell
kyleaoconnell22@gmail.com
3/20/19

python capture_trim_summary_v2.py -i .


****remove sampels with _R_ in front from all calcs
go ahead and remove the ---s to calculate 
'''
#define the user inputs 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    return parser.parse_args()

#calc averages
def get_average(array):
	ave = np.average(array)
	return round(ave,2)

#calc percentage for missing data calc
def percentage(part, whole):
	#if alignment is empty (in my case the mtDNA gblocks trimmed) then pass so it doesnt crash the script
	if whole == 0:
		pass
	else:
  		return round(100 * float(part)/float(whole),2)

#Summary routine
def summarize_alignments(in_dir):
	suffix = in_dir
	print 'in_directory = ',in_dir
	#it will run this function on both output dirs
	sum_out = suffix+'_summary.txt'
	fh_out = open(sum_out, 'a')
	fh_len = open('locus_lengths.txt','a')
	if os.path.exists(in_dir):
		os.chdir(in_dir)
		# task 1 : count files in dir to know how many loci are left after trimming
		filenum = [] #list to count number of files
		indnum = [] #count num of ind per aln
		i = 0 #count total bp of sequence captured
		j = 0 #count '-' to find % missing data across alns
		seqs = [] #count length of seqs
		namelist = []  #count unique names in files
		l = []
		namecounts = []
		
		print 'reading in all loci'
		for filetype in os.listdir('.'):
			if filetype == '.DS_store':
				os.remove(filtype)
			else:
				#count number of alignments = #loci captured across samples
				filenum.append(filetype)
					
			#use SeqIO to read in the fasta files easier
			with open(filetype, "rU") as handle:
				if filetype.endswith('.fasta'):
					#print 'processing', filetype
					seq_temp = [] #used to calc average length of each aln, rather than of each seq
					for record in SeqIO.parse(handle, "fasta"):
						if not record.id.startswith('_R_'):
							namelist.append(record.id) #write all names to namelist
										
							i = i + len(record.seq.strip('-')) #calc total bp in dataset
							seq_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
						for char in record.seq:
							if char == '-' or char == 'N': 
								j = j + 1 #calc the number of dashes, we will use this to find missing data below
					
					ave = np.average(seq_temp) #this is average of all seqs in this aln
					max_len = max(seq_temp)
					#write to seqs list above just one len for the whole locus so each locus has one entry in the final list
					seqs.append(ave)
					l.append(filetype.split('.')[0]+'\t'+str(max_len))	
		
		print 'calculating stats'			
		#iterate through the namelist		
		for name in sorted(namelist):
			#write the full counts to another list
			namestats = name + '\t' + str(namelist.count(name))
			namecounts.append(namestats)
		
		print 'calculating average len of each locus, writing to outfile'
		l = set(l)
		for locus in l:
			fh_len.write(locus + '\n')
		
		##### Write the summary output stats #####		
		#write number of loci retained to output
		print 'there are {0} files in the {1} directory'.format(len(filenum),suffix)
		fh_out.write('there are {0} files in the {1} directory'.format(len(filenum),suffix))
		fh_out.write('\n\n')
		
		#write total bp in final alns
		fh_out.write('the alignments have {0} total bps'.format(i))
		fh_out.write('\n\n')
		
		#write % missing data ('-' or 'N' chars) in final alns
		fh_out.write('the alignments have {0}% missing data across all seqs'.format(percentage(j,i)))
		fh_out.write('\n\n')
		
		#write ave len of each alignment
		fh_out.write('The average len of a locus = {0} bp'.format(get_average(seqs[1:])))
		fh_out.write('\n\n')
		fh_out.write('The range of lengths = {0} to {1} bp'.format(min(seqs),max(seqs)))
		fh_out.write('\n\n')
		
		#get unique values of counts (using set and then converting back to list for sorting)
		nameset = set(namecounts)
		#convert back to list to sort it
		namecountsuniq = sorted(list(nameset))
		
		fh_out.write('{0} samples are represented in these alignments'.format(len(namecountsuniq)))
		fh_out.write('\n\n')
		
		#write average individual per alignment
		fh_out.write('average number of samples in each alignment = {0}'.format(len(namelist)/len(filenum)*2))
		fh_out.write('\n\n')
		
		#write header for indiv coverage list
		fh_out.write('sample_name' + '\t' + 'captured_loci' + '\n')
		
		#iterate through and write the capture success per sample to the output file
		#indiv 1 \t #loci present
		for item in namecountsuniq:
			fh_out.write(item + '\n')
			#print item
			
	#if the dir doesn't exist then don't do anything
	#for example if you only trimmed with trimal not gblocks, then it will skip the gblocks summary routine
	else:
		pass
	os.chdir('..')

def main():
	#define the arguments
	args = get_args()
	summarize_alignments(args.in_dir)

if __name__ == '__main__':
    main()
