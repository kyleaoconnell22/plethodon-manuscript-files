import sys
import subprocess as sp
import os
'''
requires that bowtie2 is installed
Need to run bowtie2-build first

'''
#give path to root
fasta_directory = sys.argv[1]
#move to the base directory
os.chdir(fasta_directory)

#give path to reference genome
db = sys.argv[2]
##########################################################
#create a blast db for the final fasta sheet
'''
#assumes the following command is already run : bowtie2-build ../../ecoli.fa ecoli
'''
samples = []
for filetype in os.listdir('.'):
	if filetype.endswith('R1.fastq'):
		samples.append(filetype.split('_')[0]+'_'+filetype.split('_')[1])

#  my $call2 = system("$bowtie -x $contam -1 $reads{'1'} -2 $reads{'2'} --fast -S $contamout1 --sam-nohead --sam-nosq");


for s in samples:
	query = filetype
	outfile = query + '_'+ db + '.blastout'
	forward = s+'_R1.fastq'
	rev = s+'_R2.fastq'
	out = s+'_'+db+'.sam'
	proc_bowtie2 = sp.Popen(['bowtie2', '-x', db, '-1', forward, '-2', rev, '-S', out])
	proc_bowtie2.wait()