import numpy as np
import subprocess as sp
import os

fh_out=open('depth_summary.txt','w')

for file in os.listdir('.'):
	# Open the VCF files and grab the filename
	if file.endswith('.vcf'):
		name=file.split('.')[0]
		
		# Run VCFTOOLS to get mean depth per site		
		vcf_sub = f"vcftools --vcf {file} --site-mean-depth --out {name}"
		proc_sub = sp.call(vcf_sub,shell=True)
		
		temp=name+'.ldepth.mean'
		fh_temp=open(temp,'r')
		ave_list=[]
		for line in fh_temp:
			line=line.strip()
			if not line.startswith('CHROM'): #Skip header line
				ave_list.append(float(line.split('\t')[2]))
		
		mean=round(np.mean(ave_list),2)
		
		fh_out.write(name+'\t'+str(mean)+'\n')
				