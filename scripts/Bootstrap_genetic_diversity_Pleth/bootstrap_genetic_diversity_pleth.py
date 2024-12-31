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
    parser.add_argument("-max", "--max_ind", required=True, help="REQUIRED: Max number of individuals to randomly sample = largest dataset")
    parser.add_argument("-v", "--in_vcf", required=True, help="REQUIRED: VCF with full dataset")
    parser.add_argument("-R", "--replicates", required=True, help="REQUIRED: list of replicate types coresponding to existing keep files, see documentation")
    parser.add_argument("-s", "--species", required=True, help="REQUIRED: list of species in dataset")
    parser.add_argument("-l", "--localities", help="OPTIONAL: (T or F) use flag if you want stats split by locality as well as species and replicate")
    return parser.parse_args()

def keep_files(in_dir,in_vcf):
	args = get_args()
	os.chdir(in_dir)
	pop_list = []
	inds = []
	fh_vcf = open(in_vcf,'r') #open vcf
	for line in fh_vcf:
		#get the ind names
		if line.startswith('#CHROM'):
			names=line.split('\t')[9:] 
			for name in names:
				#make name list of all inds
				inds.append(name)
				#just name, rep
				if not args.localities:
					pop_list.append(name.split('_')[0]+'_'+name.split('_')[2])
				else: #else name , locality, rep
					pop_list.append(name.split('_')[0]+'_'+name.split('_')[1]+'_'+name.split('_')[2])

		pop_set = set(pop_list)
		for i in pop_set:
			keep = i+'_keep.txt'
			fh_temp = open(keep,'w')
			for ind in inds:
				if not args.localities:
					if ind.split('_')[0]+'_'+ind.split('_')[2]==keep.strip('_keep.txt'):
						fh_temp.write(ind+'\n')
				else:
					
					if ind.split('_')[0]+'_'+ind.split('_')[1]+'_'+ind.split('_')[2]==keep.strip('_keep.txt'):
						fh_temp.write(ind+'\n')
						
def remove_empty_keep(in_dir):
	os.chdir(in_dir)
	for filename in os.listdir('.'):
		if os.stat(filename).st_size == 0:
			print(" Removing ",filename)
			os.remove(filename)  

def subsample_loop(in_dir,boots,min_ind,max_ind,in_vcf,replicates,species):
	args = get_args()
	os.chdir(in_dir)
	sp_list = [i.strip() for i in open(species)]
	rep_list = [i.strip() for i in open(replicates)]
	loc_list = ['SG','IGG'] ###this needs to be explicity stated if -l == T
	
	if not args.localities:
		for species in sp_list:
			for rep in rep_list:
				#downsample vcf
				vcf	= in_vcf
				keep_file = species + '_' + rep+'_keep.txt'			
				outfile = species + '_'+ rep+'_sumstats_summary_combined.txt'
				fh_out = open(outfile,'w')
				fh_out.write('# Pop ID'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
				+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
				+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
				+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")

				for i in range(int(boots)): #create boots loops
					j= random.randrange(int(min_ind),int(max_ind))
					#print i,j
					out = keep_file.replace('_keep','_temp')+'.'+str(i)+'_'+str(j)
					#call vcftools to subsample randomly, max-indv randomly thins
					
					print '\n'*3
					print 'subsampling vcf to ', out
					vcf_sub = "vcftools --vcf {0} --keep {1} --max-indv {2} --recode --out {3} ".format(vcf,keep_file,j,out)
					proc_sub = sp.call(vcf_sub,shell=True)
				
					pops_in = out+'.recode.vcf'
					if os.path.exists(pops_in):
						print '\n'*3
						print 'running populations on ',pops_in
						populations_call = "populations -V {0} -O .".format(pops_in)
						proc_populations = sp.call(populations_call,shell=True)
					
						sumstats = out + '.recode' + '.p.sumstats_summary.tsv'
						last_line = open(sumstats, "r").readlines()[-1]
						print '\n'*3
						print 'writing pops info to combined ', sumstats

						#write last line to file
						fh_out.write(last_line)
				#remove temp files
				for filetype in os.listdir('.'):
					if filetype.endswith('.log') or filetype.endswith('.tsv') or 'temp' in filetype:
						os.remove(filetype)

	#split by locality
	else:
		for species in sp_list:
			for loc in loc_list:
				for rep in rep_list:
					#downsample vcf
					vcf	= in_vcf
					keep_file = species + '_' + loc + '_' + rep +'_keep.txt'			
					outfile = species +'_' + loc + '_'+ rep+'_sumstats_summary_combined.txt'
					fh_out = open(outfile,'w')
					fh_out.write('# Pop ID'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
					+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
					+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
					+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")

					for i in range(int(boots)): #create boots loops
						j= random.randrange(int(min_ind),int(max_ind))
						#print i,j
						out = keep_file.replace('_keep','_temp')+'.'+str(i)+'_'+str(j)
						#call vcftools to subsample randomly, max-indv randomly thins
					
						print '\n'*3
						print 'subsampling vcf to ', out
						vcf_sub = "vcftools --vcf {0} --keep {1} --max-indv {2} --recode --out {3} ".format(vcf,keep_file,j,out)
						proc_sub = sp.call(vcf_sub,shell=True)
				
						pops_in = out+'.recode.vcf'
						if os.path.exists(pops_in):
							print '\n'*3
							print 'running populations on ',pops_in
							populations_call = "populations -V {0} -O .".format(pops_in)
							proc_populations = sp.call(populations_call,shell=True)
					
							sumstats = out + '.recode' + '.p.sumstats_summary.tsv'
							last_line = open(sumstats, "r").readlines()[-1]
							print '\n'*3
							print 'writing pops info to combined ', sumstats

							#write last line to file
							fh_out.write(last_line)
																
					#remove temp files
					for filetype in os.listdir('.'):
						if filetype.endswith('.log') or filetype.endswith('.tsv') or 'temp' in filetype:
							os.remove(filetype)


def rename_outs(in_dir,replicates,species):
	args = get_args()
	os.chdir(in_dir)
	if not args.localities:
		fh_out = open('no_locs_sumstats_combined.txt','w')
		fh_out.write('Species'+'\t'+'Replicate'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
		+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
		+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
		+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")
	
		sp_list = [i.strip() for i in open(species)]
		rep_list = [i.strip() for i in open(replicates)]

		for species in sp_list:
			for rep in rep_list:
				temp = species + '_' + rep + '_sumstats_summary_combined.txt'
				fh_temp = open(temp,'r')
				for line in fh_temp:
					if line.startswith('#'):
						pass
					else:
						line=line.strip()
						line=line.split('\t')
						fh_out.write(species + '\t'+rep+'\t')
						for i in line[1:]:
							fh_out.write(i+'\t')
						fh_out.write('\n')


	else:
		fh_out = open('locs_sumstats_combined.txt','w')
		fh_out.write('Species'+'\t'+ 'Locality'+'\t'+'Replicate'+'\t'+'Private'+'\t'+'Sites'+'\t'+'Variant_Sites'+'\t'+"Polymorphic_Sites"+'\t'+"%Polymorphic_Loci"\
		+'\t'+"Num_Indv"+'\t'+"Var_N"+'\t'+"StdErr_N"+'\t'+"P"+'\t'+"Var_P"+'\t'+"StdErr_P"+'\t'+"Obs_Het"+'\t'+"Var_OHe"+'\t'+"StdErr_OHe"\
		+'\t'+"Obs_Hom"+'\t'+"Var_OHo"+'\t'+"StdErr_OHo"+'\t'+"Exp_Het"+'\t'+"Var_EHe"+'\t'+"StdErr_EHe"+'\t'+"Exp_Hom"+'\t'+"Var_EHo"+'\t'+"StdErr_EHo"\
		+'\t'+"Pi"+'\t'+"Var_Pi"+'\t'+"StdErr_Pi"+'\t'+"Fis"+'\t'+"Var_Fis"+'\t'+"StdErr_Fis"+"\n")
	
		sp_list = [i.strip() for i in open(species)]
		rep_list = [i.strip() for i in open(replicates)]
		loc_list = ['SG','IGG'] ###this needs to be explicity stated if -l == T
	
		for species in sp_list:
			for loc in loc_list:
				for rep in rep_list:
					temp = species +'_' + loc + '_'+ rep + '_sumstats_summary_combined.txt'
					fh_temp = open(temp,'r')
					for line in fh_temp:
						if line.startswith('#'):
							pass
						else:
							line=line.strip()
							line=line.split('\t')
							fh_out.write(species + '\t' + loc + '\t'+rep+'\t')
							for i in line[1:]:
								fh_out.write(i+'\t')
							fh_out.write('\n')

	#remove sumstats summary files
	for file in os.listdir('.'):
		if file.endswith('summary_combined.txt'):
			os.remove(file)

def main():
	#define the arguments
	args = get_args()
	#keep_files(args.in_dir,args.in_vcf)
	#remove_empty_keep(args.in_dir)
	subsample_loop(args.in_dir,args.boots,args.min_ind,args.max_ind,args.in_vcf,args.replicates,args.species)
	remove_empty_keep(args.in_dir)
	rename_outs(args.in_dir,args.replicates,args.species)

if __name__ == '__main__':
    main()	
	
	
	
	
	