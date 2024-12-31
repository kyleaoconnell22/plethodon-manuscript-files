#!/usr/bin/python

import argparse
import os
import subprocess as sp
import shutil
import numpy as np
import sys


'''

filter to each pop, greater than 3 samples
find allele freqs
open allele freq

split allele freq file into a1 and a2
k = 0
if a1 or a2 is zero then add to counter
count total alleles
total/k = perc lost alleles
fh_out write pop and percent lost
'''

def get_args():
	"""
	Get arguments from command line.
	"""
	parser = argparse.ArgumentParser(
		description="""--------------------------------------------------------------------------_
   	Script will create keep files (requires pop file and sample file)
   	then will calc allele freqs for each pop
   	and then calc perc fixed alleles at zero for each pop (indicative of allele loss)
   	
    kyleaoconnell22@gmail.com
    Jan 2021
    --------------------------------------------------------------------------_""")
	parser.add_argument("-v", "--in_vcf",required=True,help="REQUIRED: full path to vcf file with all samples")
	parser.add_argument("-s", "--samples",required=True,help="REQUIRED: list of all samples you want to include in the analysis")
	parser.add_argument("-p", "--pops",required=True,help="REQUIRED: pops you want to calc allele freqs for")

	return parser.parse_args()

fh_out = open('fixed_alleles.txt','w')
fh_out.write('population'+'\t'+'species'+'\t'+'locality'+'\t'+'replicate'+'\t'+'lost_alleles'+'\n')

def gen_keep(in_vcf,samples,pops):
	pop_list = [i.strip() for i in open(pops)]
	samps = [i.strip() for i in open(samples)]
	for p in pop_list:
		f=p+'_keep.txt'
		fh_temp = open(f,'w')
		for s in samps:
			if s.startswith(p):
				fh_temp.write(s+'\n')

def calc_freqs(in_vcf):
	for filehandle in os.listdir('.'):
		if filehandle.endswith('keep.txt'):
			call_vcftools = ("vcftools --vcf {0} --keep {1} --freq --out {2}".format(in_vcf,filehandle,'freq.temp.out'))
			proc = sp.call(call_vcftools, shell=True)
			fh_freq = open('freq.temp.out.frq','r')
			t = 0
			z = 0
			for line in fh_freq:
				line=line.strip()
				if line.startswith('CHROM'):
					pass
				else:
					a1=line.split('\t')[4].split(':')[1]
					a2=line.split('\t')[5].split(':')[1]
					t = t + 1
					if a1 == '0' or a2 == '0':
						z = z+1
			p_zero = str(round(float(z)/float(t),2))
			pops = filehandle.strip('_keep.txt')
			spec = filehandle.split('_')[0]
			loc = filehandle.split('_')[1]
			rep = filehandle.split('_')[2]
			fh_out.write(pops+'\t'+spec+'\t'+loc+'\t'+rep+'\t'+p_zero+'\n')
						
			
		
			
			
			
			
def main():
	args = get_args()
	#gen_keep(args.in_vcf,args.samples,args.pops)
	calc_freqs(args.in_vcf)
	
	
if __name__ == '__main__':
    main()