import os
import sys
import subprocess as sp
import shutil
import argparse
import random
'''
need replicate_list.txt which needs to list the replicates like this:
supernatant
pellet
formalin-fixed
RAD

And these need to correspond to existing files called:
supernatant_keep.txt
pellet_keep.txt 
etc. Which need to list all the samples for this replicate-type, and match the names in the vcffile

Species list needs to match species names in sample names in the vcf

Need to have an output that is sorted by species and datatype,
and another one that is species, locality, datatype

'''

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-b", "--boots", required=True, help="REQUIRED: Number of bootstraps.")
    parser.add_argument("-min", "--min_ind", required=True, help="REQUIRED: Min number of individuals to randomly sample = smallest dataset")
    #parser.add_argument("-max", "--max_ind", required=True, help="REQUIRED: Max number of individuals to randomly sample = largest dataset")
    parser.add_argument("-v", "--in_vcf", required=True, help="REQUIRED: VCF with full dataset")
    parser.add_argument("-r", "--replicates", required=True, help="REQUIRED: list of replicate types coresponding to existing keep files, see documentation")
    parser.add_argument("-s", "--samples", required=True, help="REQUIRED: list of samples in dataset")
    parser.add_argument("-l", "--localities", help="OPTIONAL: (T or F) use flag if you want stats split by locality as well as species and replicate")
    parser.add_argument("-p", "--pops",required=True,help="REQUIRED: pops you want to include in analysis")

    return parser.parse_args()

def gen_keep(in_vcf,samples,pops):
	#generates a keep file for each pop from pop file and writes all samples from sample file to the correct file
	#assumes names are species_locality_replicate_indiv
	#needs to be run once then can be commented out below
	pop_list = [i.strip() for i in open(pops)]
	samps = [i.strip() for i in open(samples)]
	for p in pop_list:
		f=p+'_keep.txt'
		fh_temp = open(f,'w')
		for s in samps:
			if s.startswith(p):
				fh_temp.write(s+'\n')


def subsample_loop(in_dir,boots,min_ind,in_vcf,replicates,samples,pops):
	args = get_args()
	os.chdir(in_dir) #change into base dir
	out_name=in_vcf.split('.')[0]+'sumstats_combined_all_samples.txt'
	fh_out = open(out_name,'w')
	fh_out.write('Population'+'\t'+'Species'+'\t'+'Locality'+'\t'+'Replicate'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
	+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
	+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
	+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")
	 
	pop_list = [i.strip() for i in open(pops)] #to get #pops below
	
	for n in range(int(boots)): #create boots loops
		#write subsampled keep file for each of the pops and combined into one file
		fh_keep_tp = open('keep.tmp','w') #combined keep file
		fh_pops_tp = open('pops_in.tmp','w') #combined subsampled pop map for populations
		
		for filehandle in os.listdir('.'):
			if filehandle.endswith('keep.txt'): #grab the keep files
				inds = [i.strip() for i in open(filehandle)] #get inds
				j= random.randrange(int(min_ind),len(inds)) #generate random number for each pop between min (2) and len of file (all inds)
				k = random.sample(inds,j) #randomly subsample the pop based on j
				for i in k: #write the indiv. to keep and pop map
					fh_keep_tp.write(i+'\n')
					population = i.split('_')[0]+'_'+i.split('_')[1]+'_'+i.split('_')[2] #assumes names are species_locality_replicate_indiv
					fh_pops_tp.write(i+'\t'+population+'\n')
		
		out = in_vcf.split('.')[0]+'.'+str(n)+'_temp' #vcf tmp out name
		vcf_sub = "vcftools --vcf {0} --keep {1} --recode --out {2} ".format(in_vcf,'keep.tmp',out) #subsample original vcf based on new keep file
		proc_sub = sp.call(vcf_sub,shell=True)
		
		#run populations on subsampled vcf
		pops_in = out+'.recode.vcf'
		if os.path.exists(pops_in):
			print ('\n'*3)
			print ('running populations on ',pops_in)
			populations_call = "populations -V {0} -O . -M {1}".format(pops_in,'pops_in.tmp')
			proc_populations = sp.call(populations_call,shell=True)
	
		#write populations output for each boot to the final output  file
		sumstats = out + '.recode.p.sumstats_summary.tsv'
		
		if not n == 0:
			lines = open(sumstats, "r").readlines()
			r=len(pops)+13 #to skip the header rows, but play with this and make sure you have all the pops being printed out
			k_lines = lines[r:]
			for l in k_lines:
				#print(l)
				l=l.split('\t')
				population = l[0]
				spec = population.split('_')[0]
				rep=population.split('_')[1]
				loc=population.split('_')[2]
				fh_out.write(population+'\t'+spec+'\t'+rep+'\t'+loc+'\t')
				for s in l[1:-1]:
					fh_out.write(s+'\t')
				fh_out.write(l[-1])
		
		#remove tmp files as you go to reduce mem req
		for filetype in os.listdir('.'):
			if 'temp' in filetype or 'distribs' in filetype or 'log' in filetype:
				os.remove(filetype)
		
def main():
	#define the arguments
	args = get_args()
	#gen_keep(args.in_vcf,args.samples,args.pops)
	subsample_loop(args.in_dir,args.boots,args.min_ind,args.in_vcf,args.replicates,args.samples,args.pops)

if __name__ == '__main__':
    main()	
	
	
	
	
	