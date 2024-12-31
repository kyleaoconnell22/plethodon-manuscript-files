#!/usr/bin/python

import argparse
import os
import sys
import subprocess as sp
import shutil
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np

'''
Need to calculate number/ total of each type of seq in each aln
Need to calculate missing data for each type per sequence
Need to calculate average length of aln for each data type

'''

in_dir = sys.argv[1]
os.chdir(in_dir)
fh_out = open('data_type_summary.txt','a')

#lists for finding number of samples per aln
ffs = []
bl = []
allo = []
fr = []

#list for finding # of samples per data type
ffs_names = []
bl_names = []
allo_names = []
fr_names = []

#lists for finding average leng
ffs_len = []
bl_len = []
allo_len = []
fr_len = []

#lists for finding missing data
ffs_miss = []
bl_miss = []
allo_miss = []
fr_miss = []

for filetype in os.listdir('.'):
	if filetype.endswith('.fasta'):
		
		#Plethodon data
		if filetype.startswith('set'):
			#use SeqIO to read in the fasta files easier
			with open(filetype, "rU") as handle:
				#count n and len
				ffs_temp = []
				bl_temp = []
				fr_temp = []
			
				#count missing
				ffs_m = 0
				bl_m = 0
				fr_m = 0
			
				for record in SeqIO.parse(handle, "fasta"):
					#skip empty sequences
					if len(record.seq.strip('-'))>0:
						ffs_j = 0
						bl_j = 0
						fr_j = 0
						if not record.id.startswith('_R_'): #this happens in the phased stuff
							if 'FFS' in record.id:
								#print record.id
								ffs_names.append(record.id)
								ffs_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								#print len(record.seq.strip('-'))
								#print record.seq.strip('-')
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										ffs_j = ffs_j + 1 #calc the number of dashes, we will use this to find missing data below
							elif 'Blood' in record.id:
								bl_names.append(record.id)
								bl_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										bl_j = bl_j + 1 #calc the number of dashes, we will use this to find missing data below
							elif 'fresh' in record.id or 'Fresh' in record.id:
								fr_names.append(record.id)
								fr_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										fr_j = fr_j + 1 #calc the number of dashes, we will use this to find missing data below
						ffs_m = ffs_j
						bl_m = bl_j
						fr_m = fr_j
	
						#calc averages across samples for each indiv aln
						#numbers of samples
						if len(ffs_temp) > 0:
							ffs.append(len(ffs_temp)) #write FFS sample numbers to the top list
						elif len(bl_temp)>0:
							bl.append(len(bl_temp)) #write FFS sample numbers to the top list
						elif len(fr_temp)>0:
							fr.append(len(fr_temp)) #write FFS sample numbers to the top list
	
						#length
						if len(ffs_temp) > 0:
							ave_len_ffs = float(np.average(ffs_temp))
							ffs_len.append(ave_len_ffs) #write len of aln to the top list
						
						elif len(bl_temp)>0:
							ave_len_bl = np.average(bl_temp)
							bl_len.append(ave_len_bl) #write len of aln to the top list
						elif len(fr_temp)>0:
							ave_len_fr = np.average(fr_temp)
							fr_len.append(ave_len_fr) #write len of aln to the top list
	
						#missing
						if ffs_j > 0:
							ave_miss_ffs = float(ffs_j/ave_len_ffs)
							ffs_miss.append(ave_miss_ffs)
						elif bl_j > 0:
							ave_miss_bl = float(bl_j/ave_len_bl)
							bl_miss.append(ave_miss_bl)
						elif fr_j > 0:
							ave_len_fr = np.average(fr_temp)
							ave_miss_fr = float(fr_j/ave_len_fr)
							fr_miss.append(ave_miss_fr)
							
	
		if filetype.startswith('ocus'):
			#use SeqIO to read in the fasta files easier
			with open(filetype, "rU") as handle:
				#count n and len
				ffs_temp = []
				allo_temp = []
				fr_temp = []
			
				#count missing
				ffs_m = 0
				allo_m = 0
				fr_m = 0
			
				for record in SeqIO.parse(handle, "fasta"):
					#skip empty sequences
					if len(record.seq.strip('-'))>0:
						ffs_j = 0
						allo_j = 0
						fr_j = 0
						if not record.id.startswith('_R_'): #this happens in the phased stuff
							if 'FFS' in record.id:
								#print record.id
								ffs_names.append(record.id)
								ffs_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								#print len(record.seq.strip('-'))
								#print record.seq.strip('-')
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										ffs_j = ffs_j + 1 #calc the number of dashes, we will use this to find missing data below
							elif 'allozyme' in record.id:
								allo_names.append(record.id)
								allo_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										allo_j = allo_j + 1 #calc the number of dashes, we will use this to find missing data below
							elif 'fresh' in record.id:
								fr_names.append(record.id)
								fr_temp.append(len(record.seq.strip('-'))) #write just the len of the aln, not every seq
								for char in record.seq:
									if char == '-' or char == 'N' or char == 'n': 
										fr_j = fr_j + 1 #calc the number of dashes, we will use this to find missing data below
						ffs_m = ffs_j
						allo_m = allo_j
						fr_m = fr_j
						#fr is acting up so I moved it here above the 0 loop filters below
						fr.append(len(fr_temp)) #write FFS sample numbers to the top list

						#calc averages across samples for each indiv aln
						#numbers of samples
						if len(ffs_temp) > 0:
							ffs.append(len(ffs_temp)) #write FFS sample numbers to the top list
						elif len(allo_temp)>0:
							allo.append(len(allo_temp)) #write FFS sample numbers to the top list
						#length
						if len(ffs_temp) > 0:
							ave_len_ffs = float(np.average(ffs_temp))
							ffs_len.append(ave_len_ffs) #write len of aln to the top list
						
						elif len(allo_temp)>0:
							ave_len_allo = np.average(allo_temp)
							allo_len.append(ave_len_allo) #write len of aln to the top list
						
						ave_len_fr = np.average(fr_temp)
						fr_len.append(ave_len_fr) #write len of aln to the top list
	
	
						#missing
						if ffs_j > 0:
							ave_miss_ffs = float(ffs_j/ave_len_ffs)
							ffs_miss.append(ave_miss_ffs)
						elif allo_j > 0:
							ave_miss_allo = float(allo_j/ave_len_allo)
							allo_miss.append(ave_miss_allo)
						elif fr_j > 0:
							ave_len_fr = np.average(fr_temp)
							ave_miss_fr = float(fr_j/ave_len_fr)
							fr_miss.append(ave_miss_fr)

#pleth
if len(allo) == 0:
	print 'running pleth aln', 
	#Need to calculate number/ total of each type of seq in each aln
	ave_ffs_n = round(float(np.average(ffs)/len(set(ffs_names)))*100,2)
	ave_bl_n = round(float(np.average(bl)/len(set(bl_names)))*100,2)
	ave_fr_n = round(float(np.average(fr)/len(set(fr_names)))*100,2)
	
	fh_out.write('average proportion of samples in alns for FFS = {0}, Blood = {1}, Fresh = {2}'.format(ave_ffs_n,ave_bl_n,ave_fr_n))
	fh_out.write('\n')
	
	#Need to calculate missing data for each type per sequence
	ave_ffs_miss = round(np.average(ffs_miss)*100,2)
	ave_bl_miss = round(np.average(bl_miss)*100,2)
	ave_fr_miss = round(np.average(fr_miss)*100,2)
	
	fh_out.write('average percent missing in alns for FFS = {0}, Blood = {1}, Fresh = {2}'.format(ave_ffs_miss,ave_bl_miss,ave_fr_miss))
	fh_out.write('\n')

	#Need to calculate average length of aln for each data type
	ave_ffs_len = round(np.average(ffs_len),2)
	ave_bl_len = round(np.average(bl_len),2)
	ave_fr_len = round(np.average(fr_len),2)
	
	fh_out.write('average length in alns for FFS = {0}bp, Blood = {1}bp, Fresh = {2}bp'.format(ave_ffs_len,ave_bl_len,ave_fr_len))
	fh_out.write('\n')

#gyros
else:
	print 'running gyros'
	#Need to calculate number/ total of each type of seq in each aln
	ave_ffs_n = round(float(np.average(ffs)/len(set(ffs_names)))*100,2)
	ave_allo_n = round(float(np.average(allo)/len(set(allo_names)))*100,2)
	ave_fr_n = round(float(np.average(fr)/len(set(fr_names)))*100,2)
	#print ave_fr_n
	fh_out.write('average proportion of samples in alns for FFS = {0}, Allozyme = {1}, Fresh = {2}'.format(ave_ffs_n,ave_allo_n,ave_fr_n))
	fh_out.write('\n')

	#Need to calculate missing data for each type per sequence
	ave_ffs_miss = round(np.average(ffs_miss)*100,2)
	ave_allo_miss = round(np.average(allo_miss)*100,2)
	ave_fr_miss = round(np.average(fr_miss)*100,2)
	
	fh_out.write('average % missing in alns for FFS = {0}, Allozyme = {1}, Fresh = {2}'.format(ave_ffs_miss,ave_allo_miss,ave_fr_miss))
	fh_out.write('\n')

	#Need to calculate average length of aln for each data type
	ave_ffs_len = round(np.average(ffs_len),2)
	ave_allo_len = round(np.average(allo_len),2)
	cleaned_fr_len = [x for x in fr_len if str(x) != 'nan']
	ave_fr_len = round(np.average(cleaned_fr_len),2)
	
	fh_out.write('average length in alns for FFS = {0}bp, Allozyme = {1}bp, Fresh = {2}bp'.format(ave_ffs_len,ave_allo_len,ave_fr_len))
	fh_out.write('\n')
