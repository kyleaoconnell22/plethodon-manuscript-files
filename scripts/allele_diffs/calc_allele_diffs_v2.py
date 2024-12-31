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
fh_out = open('allele_diff_stats.txt','a')

def parse_vcf(in_vcf, sample_list):
	print 'running parse vcf'
	fh_s = open(sample_list,'r')
	samples = []
	for line in fh_s:
		line=line.strip()
		samples.append(line.split('_')[0])
	samples=set(samples)
	
	for s in samples:
	
	#this first section is because I have a muscle and liver duplicate of one sample and 
	#so it gets complex for that reason, but the next chunk would apply to most datasets with 
	#just numerical naming
		if s == '525133M':
			pass
		elif s == '525133L':
			pass
		elif s == '525133':
			ffsM_al = s+'_M_ffs_al'
			ffsL_al = s+'_L_ffs_al'
			ffsM_ffsL = s+'L_M_ffs'
			ffsM_fr = s+'_M_ffs_fr'
			ffsL_fr = s+'_L_ffs_fr'
			al_fr = s + '_al_fr'
			
			fh_ffsM_al = open(ffsM_al,'a')
			fh_ffsL_al = open(ffsL_al,'a')
			fh_ML = open(ffsM_ffsL,'a')
			fh_ffsM_fr = open(ffsM_fr,'a')
			fh_ffsL_fr = open(ffsL_fr,'a')
			fh_al_fr = open(al_fr,'a')
			
			fh_ffsM_al.write(s+'M_FFS_sorted.bam'+'\n'+s+'_allozyme_sorted.bam'+'\n')
			fh_ffsL_al.write(s+'L_FFS_sorted.bam'+'\n'+s+'_allozyme_sorted.bam'+'\n')
			fh_ML.write(s+'M_FFS_sorted.bam'+'\n'+s+'L_FFS_sorted.bam'+'\n')
			fh_ffsM_fr.write(s+'M_FFS_sorted.bam'+'\n'+s+'_fresh_sorted.bam'+'\n')
			fh_ffsL_fr.write(s+'L_FFS_sorted.bam'+'\n'+s+'_fresh_sorted.bam'+'\n')
			fh_al_fr.write(s+'_allozyme_sorted.bam'+'\n'+s+'_fresh_sorted.bam'+'\n')			
			
			out_ffsM_al = s+ '_'+ 'ffsM_al'
			out_ffsM_fr = s + '_' + 'ffsM_fr'
			out_ffsL_al = s+ '_'+ 'ffsL_al'
			out_ffsL_fr = s + '_' + 'ffsL_fr'
			out_ffsM_ffsL = s + '_' + 'ML'
			out_al_fr = s + '_' + 'al_fr'

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,fh_ffsM_al,out_ffsM_al))
			proc = sp.call(call_str, shell=True)
			
			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,fh_ffsL_al,out_ffsL_al))
			proc = sp.call(call_str, shell=True)

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,fh_ML,out_ffsM_ffsL))
			proc = sp.call(call_str, shell=True)

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,ffsM_fr,out_ffsM_fr))
			proc = sp.call(call_str, shell=True)
			
			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,ffsL_fr,out_ffsL_fr))
			proc = sp.call(call_str, shell=True)

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,al_fr,out_al_fr))
			proc = sp.call(call_str, shell=True)

		else:
			#this chunk will work for just normal naming conventions
			ffs_al = s+'_ffs_al'
			ffs_fr = s+'_ffs_fr'
			al_fr = s + '_al_fr'
			
			fh_ffs_al = open(ffs_al,'a')
			fh_ffs_fr = open(ffs_fr,'a')
			fh_al_fr = open(al_fr,'a')
			
			fh_ffs_al.write(s+'_FFS_sorted.bam'+'\n'+s+'_allozyme_sorted.bam'+'\n')
			fh_ffs_fr.write(s+'_FFS_sorted.bam'+'\n'+s+'_fresh_sorted.bam'+'\n')
			fh_al_fr.write(s+'_allozyme_sorted.bam'+'\n'+s+'_fresh_sorted.bam'+'\n')			
			
			
			out_ffs_al = s+ '_'+ 'ffs_al'
			out_ffs_fr = s + '_' + 'ffs_fr'
			out_al_fr = s + '_' + 'al_fr'

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,ffs_al,out_ffs_al))
			proc = sp.call(call_str, shell=True)
			
			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,ffs_fr,out_ffs_fr))
			proc = sp.call(call_str, shell=True)

			call_str = ("vcftools --vcf {0} --keep {1} --recode --out {2} --max-missing 1".format(in_vcf,al_fr,out_al_fr))
			proc = sp.call(call_str, shell=True)

def gen_str():
	for vcf in os.listdir('.'):
		if vcf.endswith('.recode.vcf'):
			out_name = vcf.split('.')[0]
			call_str = ("plink --vcf {0} --recode structure --out {1} --double-id -allow-extra-chr".format(vcf,out_name))
			proc = sp.call(call_str, shell=True)

def find_diffs():
	len_ffs_al = []
	len_ffs_fr = []
	len_al_fr = []
	
	#lists for non het/hom diffs
	d_ffs_al = []
	d_ffs_fr = []
	d_al_fr = []
	
	#make lists for het/hom diffs
	d_ffs_al_het = []
	d_ffs_fr_het = []
	d_al_fr_het = []
	
	#make lists for total diffs
	d_ffs_al_tot = []
	d_ffs_fr_tot = []
	d_al_fr_tot = []
	

	#write header to outfile
	fh_out.write('data_set'+'\t'+'%diff_alleles_het'+'\t'+'%diff_alleles_other'+'\t'+'#diffs_total'+'\t'+'#sites'+'\n')
	for struc in os.listdir('.'): #iterate structure files
		if struc.endswith('strct_in'):
			seq1 = []
			seq2 = []
			fh_temp = open(struc,'r')
			
			#write whole sequence to two lists
			for i, line in enumerate(fh_temp):
				if i == 2:
					seq1.append(line.strip().split(' ')[1:])
				elif i == 3:
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
				
					if a1 == a3 and a2 == a4: #same either hom or het, pass
						pass
					elif a1 != a3 and a2 == a4: #het error, 1 allele is diff one shared
						j = j+1
					elif a1 == a3 and a2 != a4: #het error, other allele is diff
						j = j+1
					elif a1 != a3 and a2 != a4: #both alleles are different
						i = i+1
				#print float(i), float(j), len(s1b)

										
				diff_i = round(float(i)/float(len(s1b))*100,2)
				diff_j = round(float(j)/float(len(s1b))*100,2)
				
				diff_total = round(float(i + j)/float(len(s1b))*100,2)
				fh_out.write(struc.split('.')[0]+'\t'+str(diff_j)+'\t'+str(diff_i)+'\t'+str(diff_total)+'\t'+str(len(s1b))+'\n')
			
			
				#calculate stats by datatype
				if len(s1b) > 199:
					if 'ffs_al' in struc:
						len_ffs_al.append(len(s1b))
						d_ffs_al.append(diff_i)
						d_ffs_al_het.append(diff_j)
						d_ffs_al_tot.append(diff_total)
					elif 'ffs_fr' in struc:
						len_ffs_fr.append(len(s1b))
						d_ffs_fr.append(diff_i)
						d_ffs_fr_het.append(diff_j)
						d_ffs_fr_tot.append(diff_total)
					elif 'al_fr' in struc:
						len_al_fr.append(len(s1b))
						d_al_fr.append(diff_i)
						d_al_fr_het.append(diff_j)
						d_al_fr_tot.append(diff_total)

	#write to outfile
	fh_out.write('\n'+'\n'+'stats by data set'+'\n'+'\n')
	
	#length of shared SNPs
	ave_len_ffs_al = round(np.mean(len_ffs_al),2)
	fh_out.write('ave len ffs - al = ' + str(ave_len_ffs_al) + '['+str(min(len_ffs_al))+'\t' + str(max(len_ffs_al))+']'+'\n')
	ave_len_ffs_fr = round(np.mean(len_ffs_fr),2)
	fh_out.write('ave len ffs - fr = ' + str(ave_len_ffs_fr) + '['+str(min(len_ffs_fr))+'\t' + str(max(len_ffs_fr))+']'+'\n')
	ave_len_al_fr = round(np.mean(len_al_fr),2)
	fh_out.write('ave len al - fr = ' + str(ave_len_al_fr) + '['+str(min(len_al_fr))+'\t' + str(max(len_al_fr))+']'+ '\n')
	
	#percentage of non het diffs
	fh_out.write('\n'+'\n')
	ave_d_ffs_al = round(np.mean(d_ffs_al),2)
	fh_out.write('ave non het diff ffs - al = '+ str(ave_d_ffs_al) +'['+str(min(d_ffs_al))+'\t' + str(max(d_ffs_al))+']'+ '\n')
	ave_d_ffs_fr = round(np.mean(d_ffs_fr),2)
	fh_out.write('ave non het diff ffs - fr = '+ str(ave_d_ffs_fr) +'['+str(min(d_ffs_fr))+'\t' + str(max(d_ffs_fr))+']'+ '\n')
	ave_d_al_fr = round(np.mean(d_al_fr ),2)
	fh_out.write('ave non het diff al - fr = '+ str(ave_d_al_fr) +'['+str(min(d_al_fr))+'\t' + str(max(d_al_fr))+']'+ '\n')
	
	#percentage of het diffs
	fh_out.write('\n'+'\n')
	ave_d_ffs_al_het = round(np.mean(d_ffs_al_het),2)
	fh_out.write('ave het diff ffs - al = '+ str(ave_d_ffs_al_het) +'['+str(min(d_ffs_al_het))+'\t' + str(max(d_ffs_al_het))+']'+ '\n')
	ave_d_ffs_fr_het = round(np.mean(d_ffs_fr_het),2)
	fh_out.write('ave het diff ffs - fr = '+ str(ave_d_ffs_fr_het) +'['+str(min(d_ffs_fr_het))+'\t' + str(max(d_ffs_fr_het))+']'+ '\n')
	ave_d_al_fr_het = round(np.mean(d_al_fr_het),2)
	fh_out.write('ave het diff al - fr = '+ str(ave_d_al_fr_het) +'['+str(min(d_al_fr_het))+'\t' + str(max(d_al_fr_het))+']'+ '\n')
	
	#percentage of total diffs
	fh_out.write('\n'+'\n')
	ave_d_ffs_al_tot = round(np.mean(d_ffs_al_tot),2)
	fh_out.write('ave total diff ffs - al = '+ str(ave_d_ffs_al_tot) +'['+str(min(d_ffs_al_tot))+'\t' + str(max(d_ffs_al_tot))+']'+ '\n')
	ave_d_ffs_fr_tot = round(np.mean(d_ffs_fr_tot),2)
	fh_out.write('ave total diff ffs - fr = '+ str(ave_d_ffs_fr_tot) +'['+str(min(d_ffs_fr_tot))+'\t' + str(max(d_ffs_fr_tot))+']'+ '\n')
	ave_d_al_fr_tot = round(np.mean(d_al_fr_tot),2)
	fh_out.write('ave total diff al - fr = '+ str(ave_d_al_fr_tot) +'['+str(min(d_al_fr_tot))+'\t' + str(max(d_al_fr_tot))+']'+ '\n')

def clean_dir():
	pass
	#remove recode.vcf
	#remove logs
	#remove structures
					
def main():
	args = get_args()
	#parse_vcf(args.in_vcf, args.samples)
	#gen_str()
	find_diffs()
if __name__ == '__main__':
    main()