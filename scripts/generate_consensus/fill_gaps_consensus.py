import os
import sys
from Bio import SeqIO
'''
python fill_gaps_consensus.py
Needs the outputs from create_sp_consensus.py in the dir. 
Need to hardcode reference species to use for filling in missing loci
Script will find what loci are present, then grab sequence from the selected reference and create 
a chimeric pseudoreference for each species
Kyle O'Connell 
June 2020
'''
loci= []
#generate locus list
reference = 'yonahlossee'
yon_in = reference + '_consensus.fasta'
fh = open(yon_in,'r')
yon_seqs = []

for line in fh:
	line = line.strip()
	if line.startswith('>'):
		name = line.strip('>')
		loci.append(name)

yon = SeqIO.parse(yon_in,'fasta')
for record in yon:
	yon_seqs.append(record.id+' '+str(record.seq))
for fasta in os.listdir('.'):
	if 'expanded' in fasta:
		os.remove(fasta)
	elif not fasta.startswith('yonahlossee'):
		species = fasta.split('_')[0]
		present = []
		fh = open(fasta,'r')
		new_out = fasta.split('.')[0]+'_expanded.fasta'
		fh_sp = open(new_out,'a')
		for line in fh:
			line=line.strip()
			fh_sp.write(line+'\n')
			if line.startswith('>'):
				present.append(line.strip('>'))
		missing = set(loci)-set(present)
		#print species, '\n', missing
		fh.close()
		#print species, '\n', missing
		
		for item in yon_seqs:
			if item.split(' ')[0] in missing:
				fh_sp.write('>'+item.split(' ')[0]+'\n'+item.split(' ')[1]+'\n')
			