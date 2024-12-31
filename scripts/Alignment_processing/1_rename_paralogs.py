import os
import sys

#need to give path to 'reference_fasta_header_info.txt'
fh_in = open(sys.argv[1],'r')
#need to give path to paralog list
fh_pa = open(sys.argv[2],'r')
fh_out = open('paralog_list_renamed.txt','a')
header = []
paralog = []
for line in fh_in:
	line=line.strip()
	header.append(line)
	
for line in fh_pa:
	line=line.strip()
	paralog.append(line)

for locus in paralogs:
	for i in header:
		if locus == i.split('\t')[0]:
			newname = i.split('\t')[1].split('|')[0].replace('/','-')
			fh_out.write(newname+'\n')


