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
ref_in = 'Salamander_seqs.fasta'
mapped = 'pleth_ref_seqs.fasta'

fh = open(ref_in,'r')
ref_seqs = []

for line in fh:
	line = line.strip()
	if line.startswith('>'):
		if '|' in line:
			name = line.strip('>').split('|')[0]
			loci.append(name)
		else:
			name=line.strip('>')
			loci.append(name)
#print loci
ref = SeqIO.parse(ref_in,'fasta')
for record in ref:
	ref_seqs.append(record.id+' '+str(record.seq))
for fasta in os.listdir('.'):
	if 'expanded' in fasta:
		os.remove(fasta)
	elif fasta == mapped:
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
		print 'number of missing loci in contig file', len(missing)
		#print species, '\n', missing
		fh.close()
		
		for item in ref_seqs:
			if '|' in item:
				if item.split(' ')[0].split('|')[0] in missing:
					fh_sp.write('>'+item.split(' ')[0].split('|')[0]+'\n'+item.split(' ')[1]+'\n')
			elif not '|' in item:
				if item.split(' ')[0] in missing:
					fh_sp.write('>'+item.split(' ')[0]+'\n'+item.split(' ')[1]+'\n')
			
