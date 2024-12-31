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

for line in fh_conf:
	line = line.strip()
	if line.startswith('['):
		pass
	elif ':_' in line:
		samples.append(line)
	elif line.startswith('i7'):
		bar7.append(line)
	elif line.startswith('i5'):
		bar5.append(line)
for a,b,c in zip(samples,bar7,bar5):
	input_list.append(a+' ' + b + ' ' + c)

###create way of grabbing the dup samples so I can come back to them
for sample in input_list:
	s_name = sample.split(' ')[0].split(':')[0]
	if s_name not in extra_list:
		extra_list.append(s_name)
		input_list_filt.append(sample)
	else:
		dup_list.append(sample)
fh_dup = open('dup_samples.txt','a')
for sample in dup_list:
	fh_dup.write(sample+'\n')
		
###then back to the old flow
if not os.path.exists('./2_cleaned_trimmed_reads'):
	os.mkdir('2_cleaned_trimmed_reads')

for sample in input_list_filt:
	s_name = sample.split(' ')[0].split(':')[0]
	#i7_new = i7.replace('*', sample.split(' ')[1]) #for testing
	#i5_new = i5.replace('*', sample.split(' ')[2]) #for testing
	#fh_out.write(sample.split(' ')[0] + ','+i7_new+','+ i5_new+'\n') #for testing
	b7 = sample.split(' ')[1].split(':')[1]
	i7_new = i7.replace('*', b7)
	i5_new = i5.replace('*', sample.split(' ')[2].split(':')[1])
	
	forward = '1_raw_fastq/{0}'.format(s_name+'_R1.fastq')
	reverse = '1_raw_fastq/{0}'.format(s_name+'_R2.fastq')
	outdir = './2_cleaned_trimmed_reads/{0}'.format(s_name)

	if os.path.exists(forward):
		if not os.path.exists('./2_cleaned_trimmed_reads/{0}'.format(s_name)): #create abc dir if it does not exist
			outdir = './2_cleaned_trimmed_reads/{0}'.format(s_name)
			os.mkdir(outdir)
		#then call trim galore with sample and -a and -a2 
		call_string = "trim_galore --phred33 --dont_gzip \
			-o {0} --paired -a {1} \
			-a2 {2} --retain_unpaired \
			{3} {4}"\
			.format(outdir,i7_new,i5_new,forward,reverse)
		proc = sp.call(call_string, shell=True)


dup_out = './duplicated_dirs'
if not os.path.exists(dup_out):
	os.mkdir(dup_out)
		
#go back and retrim the dup samples with the second barcode from the dup list
#iterate through the list
#go into each outdir from the first pass as the inputs
for sample in dup_list:
	s_name = sample.split(' ')[0].split(':')[0]
	#i7_new = i7.replace('*', sample.split(' ')[1]) #for testing
	#i5_new = i5.replace('*', sample.split(' ')[2]) #for testing
	#fh_out.write(sample.split(' ')[0] + ','+i7_new+','+ i5_new+'\n') #for testing
	i7_new = i7.replace('*', sample.split(' ')[1].split(':')[1])
	i5_new = i5.replace('*', sample.split(' ')[2].split(':')[1])
	
	forward = './2_cleaned_trimmed_reads/{0}/{1}'.format(s_name,s_name+'_R1_val_1.fq')
	reverse = './2_cleaned_trimmed_reads/{0}/{1}'.format(s_name,s_name+'_R2_val_2.fq')
	outdir = './2_cleaned_trimmed_reads/{0}'.format(s_name+'_dup')
	old_out = './2_cleaned_trimmed_reads/{0}'.format(s_name)
	
	if os.path.exists(forward):
		if not os.path.exists('./2_cleaned_trimmed_reads/{0}'.format(s_name+'_dup')): #create abc dir if it does not exist
			outdir = './2_cleaned_trimmed_reads/{0}'.format(s_name+'_dup')
			os.mkdir(outdir)

		#then call trim galore with sample and -a and -a2 
		call_string = "trim_galore --phred33 --dont_gzip \
			-o {0} --paired -a {1} \
			-a2 {2} --retain_unpaired \
			{3} {4}"\
			.format(outdir,i7_new,i5_new,forward,reverse)
		proc = sp.call(call_string, shell=True)

		shutil.move(old_out,dup_out)
	
		shutil.move(outdir,old_out)
	


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
