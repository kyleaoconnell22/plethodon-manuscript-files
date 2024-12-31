import os
from Bio import SeqIO

fh_out = open('contam.txt','a')
fh_out.write('sample'+'\t'+'reads'+'contam'+'prop'+'\n')
samples = []
for filetype in os.listdir('.'):
	if filetype.endswith('R1.fastq'):
		samples.append(filetype.split('_')[0]+'_'+filetype.split('_')[1])

for s in samples:
	#find len of fastq with bioseq
	in_fastq = s+'_R1.fastq'
	fastq_dict = SeqIO.index(in_fastq, "fastq")	
	reads_f = len(fastq_dict)
	print reads_f

	#find len of sam with list
	sam_r = []
	fh_sam = open(s+'_ecoli.sam','r')
	for line in fh_sam:
		if not line.startswith('@'):
			sam_r.append(line.split(':')[0])
	reads_s = len(sam_r)
	print reads_s
	
	#find prop
	prop = float(reads_s/reads_f)
	
	fh_out.write(s+'\t'+str(reads_f) + '\t' + str(reads_s) + str(prop)+'\n')
