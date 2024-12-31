import os
import sys

locus_list=[]
fh_out = open('potential_paralogs.txt','a')
i = 0
for dir in os.listdir('.'):
	if not dir.endswith('.fasta') and not dir.endswith('.txt'):
		#os.chdir(dir)
		infile = dir+'/info_paralogous_loci.txt'
		fh = open(infile,'r')
		i = i + 1
		for line in fh:
			line=line.strip()
			#grab loci with at least 4 contigs matching
			if len(line.split('\t')) > 3:
				locus_list.append(line.split('\t')[0])
		
locus_set = set(locus_list)

for locus in locus_set:
	if locus_list.count(locus) > .75*i:
		fh_out.write(locus+'\n')
	
			