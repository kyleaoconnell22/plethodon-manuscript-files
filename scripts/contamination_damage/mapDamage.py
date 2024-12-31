import os
import sys
import shutil
import argparse

'''
in dir will have sub directories, which will then each have subdir of samples
Right now only the mapdamage section works, but if you want to call a consensus before phasing that could be done
Make sure to give full paths to all inputs
-i = in dir, needs to be full path
-s species list, list of species. If you don't have multiple species give the name of the dir that has the individuals inside
-r not actually used, but could be edited to copy the pdf outputs from all runs to this dir
-R ref_dir: this needs to be a dir that has the reference used for mapping, and needs to be named species_consensus.fasta
	where species = specie in the species_list.txt
-c min coverage is not used either right now, but would be needed for calling a consensus of the scaled bam file
#example command
mapDamage -i gyro/525249_FFS_remapped/525249_FFS_no_dupls_sorted.bam -r ../fasta_ref/target_loci_gyro.fasta -n 500000 --rescale --merge-reference-sequencesmapDamage -i gyro/525249_FFS_remapped/525249_FFS_no_dupls_sorted.bam -r ../fasta_ref/target_loci_gyro.fasta -n 500000 --rescale --merge-reference-sequences
Kyle O'Connell
July 2020
'''
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains species subdirs")
    parser.add_argument("-s", "--sp_list", required=True, help="REQUIRED: The full path to a species list (or dataset list)")
    parser.add_argument("-r", "--results_dir", required=True, help="REQUIRED: The full path to the results dir")
    parser.add_argument("-R", "--ref_dir", required=True, help="REQUIRED:The full path to the dir with the reference fasta")
    #parser.add_argument("-c", "--min_cov", required=True, help="REQUIRED: Specify the min coverage to call from the bam file")
    return parser.parse_args()


def run_mapDamage(in_dir,sp_list,ref_dir,results_dir):
	os.chdir(in_dir)
	sp = [i.strip() for i in open(sp_list)]
	for species in sp:
		os.chdir(species) #change into each species dir
		for dir in os.listdir('.'): #iterate through the ind dirs
			if 'FFS' in dir: #only grab the FFS dirs
				os.chdir(dir) #enter each FFS dir
				bam = dir.strip('remapped')+'no_dupls_sorted.bam' #name of bam file
				fix = dir.strip('remapped')+'no_dupls_sorted_unscaled.bam'
				if os.path.exists(fix):
					os.rename(fix,bam)
				ref_path = ref_dir+'/'+species+'_consensus.fasta' #path to the reference
				result_path = 'results_'+bam.split('.')[0] #path to mapDamage results dir
				#mapDamage command
				print 'running mapDamage2 on ', bam
				md_cmd = 'mapDamage -i %s -r %s -n 500000 --rescale --merge-reference-sequences' %(bam,ref_path)
				os.system(md_cmd) #run it
				out_file = bam.split('.')[0]+'.rescaled.bam' #path to mapDamage rescaled bam
				unscaled = bam.split('.')[0]+'_unscaled.bam' #path to new name for unscaled bam
				os.rename(bam,unscaled) #rename secapr bam to unscaled
				
				result = result_path+'/'+out_file
				print 'moving', result
				shutil.move(result,'.')
				
				os.rename(out_file,bam) #rename scaled bam to secapr name for phasing
				joined_fasta = dir.strip('remapped')+'bam_consensus.fasta' #path to original fasta file
				if os.path.exists(joined_fasta):
					os.remove(joined_fasta) #remove fasta file
				os.chdir('..')
	#os.chdir('..')
def main():
	#define the arguments
	args = get_args()
	run_mapDamage(args.in_dir,args.sp_list,args.ref_dir,args.results_dir)

if __name__ == '__main__':
    main()
    