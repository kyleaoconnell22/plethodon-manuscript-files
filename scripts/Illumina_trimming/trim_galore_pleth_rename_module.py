import os
import subprocess as sp
import shutil  

'''
if the file has multiple barcodes, I could trim each of the original abc files sep
or, I could trim the cat file once, then go into the output and trim it again, 
then have to somehow rename these files


For trimming: 
trim AGATCGGAAGAGC
The naming goes i7-i5
Need to write a script to generate the trimmomatic file that will match the sample names from the file list to the samples1 and 2 files and then list the i7 and i5 like this:
i7-AES1158:ACACGGTT
i5-AES1158:ACACGGTT

module load bioinformatics/trim_galore/0.6.4	

i7 is forward
i5 is reverse

trim_galore --phred33 --illumina --dont_gzip -o ../2_cleaned_trimmed_reads/ --paired 525132_allozyme_R1.fastq 525132_allozyme_R2.fastq --retain_unpaired
Need to write script that uses the contig file and grabs the sample name the i7 and i5 barcodes
and then runs trim galore with the full barcode seq from the top plus the i5 and i7

then run the command as follows 
trim_galore --phred33 --dont_gzip -o ../2_cleaned_trimmed_reads/ --paired -a fulli7 -a2 fulli5 525132_allozyme_R1.fastq 525132_allozyme_R2.fastq --retain_unpaired
Then cd into outdirs and rename outputs

#trimgalore outs examples
525132_allozyme_R1_val_1.fq
525132_allozyme_R2_val_2.fq
525132_allozyme_R1_unpaired_1.fq
525132_allozyme_R2_unpaired_1.fq

#trimmomatic outs examples
525134_FFS_clean-READ1.fastq
525134_FFS_clean-READ2.fastq
525134_FFS_clean-READ1-single.fastq
525134_FFS_clean-READ2-single.fastq

input list goes name, i7, i5
yonahlossee_IGG_FFS_461352:_ 
i7-yonahlossee_IGG_FFS_461352:TGCATACA 
i5-yonahlossee_IGG_FFS_461352:TTAACACA


********make sure job script calls:
bioinformatics/trim_galore/0.6.4
And I move an updated version of Config_cleaning without the i5 and i7 universals at the top!

'''

i7 = 'CCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG'
i5 = 'AATGATACGGCGACCACCGAGATCTACAC*TCGTCGGCAGCGTC'

config = 'trimmomatic_config.txt'
fh_conf = open(config,'r')
fh_out = open('testing.txt','a')
samples = []
bar7 = []
bar5 = []
input_list = []
extra_list = []
input_list_filt = []
dup_list = []



#then after all trim galore is done, go into each outdir and rename them to match the trimmomatic outfiles
os.chdir('./2_cleaned_trimmed_reads')

for dir in os.listdir('.'):
	os.chdir(dir)
	sample = dir
	dup = sample+'_R1_val_1_val_1.fq'
	non_dup = sample+'_R1_val_1.fq'
	if os.path.exists(dup):	
		os.rename(sample+'_R1_val_1_val_1.fq',sample+'_clean-READ1.fastq') #paired R1
		os.rename(sample+'_R2_val_2_val_2.fq',sample+'_clean-READ2.fastq') #paired R2
		os.rename(sample+'_R1_val_1_unpaired_1.fq',sample+'_clean-READ1-single.fastq') #single R1
		os.rename(sample+'_R2_val_2_unpaired_2.fq',sample+'_clean-READ2-single.fastq') #single R1
		os.chdir('..')
	elif os.path.exists(non_dup):
		os.rename(sample+'_R1_val_1.fq',sample+'_clean-READ1.fastq') #paired R1
		os.rename(sample+'_R2_val_2.fq',sample+'_clean-READ2.fastq') #paired R2
		os.rename(sample+'_R1_unpaired_1.fq',sample+'_clean-READ1-single.fastq') #single R1
		os.rename(sample+'_R2_unpaired_2.fq',sample+'_clean-READ2-single.fastq') #single R1
		os.chdir('..')
	else:
		os.chdir('..')
