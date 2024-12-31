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
    Calc allelic heterozygosity - iterates through indiv in the vcf and then calcs the proportion of alleles that are het and those that are hom
    but the perc is based on present loci, not including missing loci
    
    kyleaoconnell22@gmail.com
    Jan 2021
    --------------------------------------------------------------------------_""")
    
	parser.add_argument("-v", "--in_vcf",
						required=True,
						help="REQUIRED: full path to vcf file with all samples")
	parser.add_argument("-n", "--num_samps",
						required=True,
						help="REQUIRED: number of samples")

	return parser.parse_args()
import string
ALPHA = string.ascii_letters

fh_out = open('perc_het_stats.txt','w')
fh_out.write('sample'+'\t'+'Species'+'\t'+'Locality'+'\t'+'Replicate'+'\t'+'total_snps'+'\t'+'num_miss'+'\t'+'perc_het'+'\t'+'perc_hom'+'\n')
def parse_vcf(in_vcf,num_samps):
	fh=open(in_vcf)
	seqs = []
	for line in fh:
		line = line.strip()
		if line.startswith('##'):
			pass
		else:
			seqs.append(line)
	
	for i in range(9,int(num_samps)+9):
		het = 0
		hom = 0 #major allele
		miss = 0
		for item in seqs:
			val = item.split('\t')[i]
			if val.startswith(tuple(ALPHA)):
				fh_out.write(val+'\t'+val.split('_')[0]+'\t'+val.split('_')[1]+'\t'+val.split('_')[2]+'\t')
				print val
			else:
				locus = val.split(':')[0]
				#fh_out.write('\n')
				if locus == './.':
					miss = miss + 1
				elif locus == '0/1' or locus == '1/0':
					het = het+1
				elif locus == '0/0' or locus == '1/1':
					hom = hom+1
					
				total = int(het + hom + miss)
		print total, miss, het, hom		
		p_het = round(float(int(het)/float((total-miss))),2)*100
		p_hom = round(float(int(hom)/float((total-miss))),2)*100
		#print p_het, p_hom, p_het+p_hom
		fh_out.write(str(total)+'\t'+str(miss)+'\t'+str(p_het)+'\t'+str(p_hom)+'\n')
			
			
			
			
def main():
	args = get_args()
	parse_vcf(args.in_vcf,args.num_samps)
	
	
if __name__ == '__main__':
    main()