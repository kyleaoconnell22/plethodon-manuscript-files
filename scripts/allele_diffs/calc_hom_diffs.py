#!/usr/bin/python

import argparse
import os
import subprocess as sp
import shutil
import numpy as np
import sys

def get_args():
	"""
	Get arguments from command line.
	"""
	parser = argparse.ArgumentParser(
		description="""--------------------------------------------------------------------------_
    Calc allele differences - will create new vcf for all samples in samples.txt, then will identify 
    triads of sample replicates. It will filter each triad for only sites present between those two 
    samples, then calculate the number of differences between each site and divides by the total number 
    of sites then writes the individual diffs to a file, as well as calcs the ave diffs for triad types.
    Also make sure that the names match in the vcf and sample file...
    The script is not that friendly or annotated well, so email me if you want to use it but have issues: 
    ****for some reason it never works the first time, so run it twice and it should work###
    kyleaoconnell22@gmail.com
    Aug. 2020
    --------------------------------------------------------------------------_""")
    
	parser.add_argument("-v", "--in_vcf",
						required=True,
						help="REQUIRED: full path to vcf file with all samples")
	parser.add_argument("-s", "--samples",
						required=True,
						help="REQUIRED: The full path to an list of all samples")
	return parser.parse_args()

fh_out = open('hom_diff_stats.txt','w')

def gen_keeps(in_vcf, sample_list):
	print 'running parse vcf'
	#for species in species list:
	samples = [i.strip() for i in open(sample_list)]
	samples=set(samples)
	
	#empty lists for sample names
	ffs = []
	al = []
	fr = []
	rad = []
	
	#write samples to the correct replicate list
	for s in samples:
		#sort samples to correct list
		if 'FFS' in s:
			ffs.append(s)
		elif 'allozyme' in s:
			al.append(s)
		elif 'fresh' in s:
			fr.append(s)
		elif 'rad' in s:
			rad.append(s)
	
	#make a list of lists	
	reps = []
	reps.append(ffs)
	reps.append(al)
	reps.append(fr)
	reps.append(rad)
	#then make loop of for i in range(0,len(list)+1):

	for rep in reps: #iterate through the three lists
		#print rep
		for i in range(0,len(rep)-1): #iterate through samples in each list
			for j in range(0,len(rep)-1):
				#name this by replicate type
				temp = rep[i].split('_')[1]+'_'+str(i) + '_' + str(j+1) + '_keep'
				fh_temp = open(temp,'w')
				fh_temp.write(rep[i]+'\n'+rep[j+1])
				#print rep[i], rep[j+1], i, j+1

	#remove duplicate comparisons and single sample lists
	pairs = []
	for filetype in os.listdir('.'):
		if filetype.endswith('keep'):
			#print filetype
			repl = filetype.split('_')[0]
			j = filetype.split('_')[1]
			k = filetype.split('_')[2]
			if not repl+'_'+j+'_'+k in pairs and not repl+'_'+k+'_'+j in pairs and not j==k:
				pairs.append(repl+'_'+j+'_'+k)
				print 'keeping', filetype
			else:
				os.remove(filetype)
				
def parse_vcf(in_vcf):
	for filetype in os.listdir('.'):
		if filetype.endswith('keep'):
			out = filetype.strip('_keep')+'_temp'
			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,filetype,out))
			proc = sp.call(call_str, shell=True)
			os.remove(filetype) #remove keep files


def gen_str():
	for filetype in os.listdir('.'):
		if filetype.endswith('_temp.recode.vcf'):
			out_name = filetype.split('.')[0]
			call_str = ("plink --vcf {0} --recode structure --out {1} --double-id -allow-extra-chr".format(filetype,out_name))
			proc = sp.call(call_str, shell=True)
			os.remove(filetype) #remove temp vcfs

def find_diffs():
	#lists of lens
	len_ffs = []
	len_al = []
	len_fr = []
	len_rad = []
	
	#lists for non het/hom diffs
	d_ffs = []
	d_al = []
	d_fr = []
	d_rad = []

	#make lists for het/hom diffs
	d_ffs_het = []
	d_al_het = []
	d_fr_het = []
	d_rad_het = []
	
	#write header to outfile
	fh_out.write('data_set'+'\t'+'diff_alleles_het'+'\t'+'diff_alleles_other'+'\t'+'#diffs_total'+'\t'+'het1'+'\t'+'het2'+'\t'+'#sites'+'\n')
	for struc in os.listdir('.'): #iterate structure files
		if struc.endswith('strct_in'):
			seq1 = []
			seq2 = []
			fh_temp = open(struc,'r')
			#print struc
			#write whole sequence to two lists
			for i, line in enumerate(fh_temp):
				if i == 2:
					#print allele 1 replicate
					#print line.split(' ')[0]
					seq1.append(line.strip().split(' ')[1:])
				elif i == 3:
					#print allele 2 replicate
					#print line.split(' ')[0]
					seq2.append(line.strip().split(' ')[1:])
			#write bp to list in order, but do them two by two
			s1 = []
			s2 = []
			s1b = []
			s2b = []
			for item in seq1:
				for bp in item:
					s1.append(bp)
			
			for item in seq2:
				for bp in item:
					s2.append(bp)
			#print len(s1), len(s2)
			#filter out failed samples
			if len(s1) > 199 and len(s2) > 199:				
				i = 0 #differences not related to het/hom
				j = 0 #differences related to het/hom
				k = 0 #rep 1 is het, rep2 is hom
				l = 0 #rep 1 is hom, rep2 is het
				#print for testing
				#print struc
				#print len(s1)+1, len(s2)+1 #sanity check for value of y (can't excede the len of seq)

				#create sliding window iters
				y = 1				
				for p in range(0,(len(s1)+1)): #iterate over all values in list (seq) 
					if y<len(s1): #if still smaller than len of seq
						#print y-1,y
						s1b.append(s1[y-1]+','+s1[y])
						s2b.append(s2[y-1]+','+s2[y])
						y = y + 2
			
				#now compare alleles from seq1 to seq2
				for a,b in zip(s1b,s2b):
					#print a, b
					a1 = a.split(',')[0]
					a2 = a.split(',')[1]
					a3 = b.split(',')[0]
					a4 = b.split(',')[1]
					
			
					#does not assume allele order so calcs both possible order for both reps
					#both hom
					if a1 == a2 and a3 == a4: #same either hom or het, pass
						if a1==a3:
							pass #same hom pass
						elif a1!=a3:
							i = i + 1
					#both het
					elif a1 != a2 and a3 != a4: #both het than also pass if allele order is not the same
						pass
					#one is het, one is hom
					elif a1 != a2 and a3 == a4: #if rep 1 is het but rep 2 is homo
						j = j + 1
						k = k + 1
					elif a1 == a2 and a3 != a4: #if rep1 is hom but rep2 is het
						j = j + 1
						l = l + 1
				#hom/het diffs	
				diff_i = round(float(i)/float(len(s1b))*100,2)
				diff_j = round(float(j)/float(len(s1b))*100,2)
			
				#calculate stats by datatype
				if len(s1b) > 199:
					fh_out.write(struc.split('.')[0]+'\t'+str(diff_j)+'\t'+str(diff_i)+'\t'+str(len(s1b))+'\n')
					if 'FFS' in struc:
						len_ffs.append(len(s1b))
						d_ffs.append(diff_i)
						d_ffs_het.append(diff_j)
					elif 'allozyme' in struc:
						len_al.append(len(s1b))
						d_al.append(diff_i)
						d_al_het.append(diff_j)
					elif 'fresh' in struc:
						len_fr.append(len(s1b))
						d_fr.append(diff_i)
						d_fr_het.append(diff_j)
					elif 'rad' in struc:
						len_rad.append(len(s1b))
						d_rad.append(diff_i)
						d_rad_het.append(diff_j)

	#write to outfile
	fh_out.write('\n'+'\n'+'stats by data set'+'\n'+'\n')
	
	#print sample_sizes
	fh_out.write('samples in ffs = ' + str(len(len_ffs))+'\n')
	fh_out.write('samples in al = ' + str(len(len_al))+'\n')
	fh_out.write('samples in fr = ' + str(len(len_fr))+'\n')
	fh_out.write('samples in rad = ' + str(len(len_rad))+'\n')

	#length of shared SNPs
	fh_out.write('\n'+'\n')
	ave_len_ffs = round(np.mean(len_ffs),2)
	fh_out.write('ave len ffs = ' + str(ave_len_ffs) + '['+str(min(len_ffs))+'\t' + str(max(len_ffs))+']'+ '  sd= ' + str(round(np.std(len_ffs),2)) + '\n')
	ave_len_al = round(np.mean(len_al),2)
	fh_out.write('ave len al = ' + str(ave_len_al) + '['+str(min(len_al))+'\t' + str(max(len_al))+']'+ '  sd= ' + str(round(np.std(len_al),2)) + '\n')
	ave_len_fr = round(np.mean(len_fr),2)
	fh_out.write('ave len fr = ' + str(ave_len_fr) + '['+str(min(len_fr))+'\t' + str(max(len_fr))+']'+ '  sd= ' + str(round(np.std(len_fr),2)) + '\n')
	ave_len_rad = round(np.mean(len_rad),2)
	fh_out.write('ave len rad = ' + str(ave_len_rad) + '['+str(min(len_rad))+'\t' + str(max(len_rad))+']'+ '  sd= ' + str(round(np.std(len_rad),2)) + '\n')

	#percentage of non het diffs
	fh_out.write('\n'+'\n')
	ave_d_ffs = round(np.mean(d_ffs),2)
	fh_out.write('ave non het diff ffs  = '+ str(ave_d_ffs) +'['+str(min(d_ffs))+'\t' + str(max(d_ffs))+']'+ '  sd= ' + str(round(np.std(d_ffs),2)) + '\n')
	ave_d_al = round(np.mean(d_al),2)
	fh_out.write('ave non het diff al = '+ str(ave_d_al) +'['+str(min(d_al))+'\t' + str(max(d_al))+']'+ '  sd= ' + str(round(np.std(d_al),2)) + '\n')
	ave_d_fr = round(np.mean(d_fr ),2)
	fh_out.write('ave non het diff fr = '+ str(ave_d_fr) +'['+str(min(d_fr))+'\t' + str(max(d_fr))+']'+ '  sd= ' + str(round(np.std(d_fr),2)) + '\n')
	fh_out.write('\n'+'\n')
	ave_d_rad = round(np.mean(d_rad),2)
	fh_out.write('ave non het diff rad = '+ str(ave_d_rad) +'['+str(min(d_rad))+'\t' + str(max(d_rad))+']'+ '  sd= ' + str(round(np.std(d_rad),2)) + '\n')

	#percentage of het diffs
	fh_out.write('\n'+'\n')
	ave_d_ffs_het = round(np.mean(d_ffs_het),2)
	fh_out.write('ave het diff ffs = '+ str(ave_d_ffs_het) +'['+str(min(d_ffs_het))+'\t' + str(max(d_ffs_het))+']'+ '  sd= ' + str(round(np.std(d_ffs_het),2)) + '\n')
	ave_d_al_het = round(np.mean(d_al_het),2)
	fh_out.write('ave het diff al = '+ str(ave_d_al_het) +'['+str(min(d_al_het))+'\t' + str(max(d_al_het))+']'+ '  sd= ' + str(round(np.std(d_al_het),2)) + '\n')
	ave_d_fr_het = round(np.mean(d_fr_het),2)
	fh_out.write('ave het diff  fr = '+ str(ave_d_fr_het) +'['+str(min(d_fr_het))+'\t' + str(max(d_fr_het))+']'+ '  sd= ' + str(round(np.std(d_fr_het),2)) + '\n')
	fh_out.write('\n'+'\n')
	ave_d_rad_het = round(np.mean(d_rad_het),2)
	fh_out.write('ave het diff  rad = '+ str(ave_d_rad_het) +'['+str(min(d_rad_het))+'\t' + str(max(d_rad_het))+']'+ '  sd= ' + str(round(np.std(d_rad_het),2)) + '\n')


def clean_dir():
	for filetype in os.listdir('.'):
		if 'temp' in filetype:
			os.remove(filetype)
def main():
	args = get_args()
	gen_keeps(args.in_vcf, args.samples)
	parse_vcf(args.in_vcf)
	gen_str()
	find_diffs()
	clean_dir()

if __name__ == '__main__':
    main()