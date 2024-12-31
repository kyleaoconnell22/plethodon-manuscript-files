'''
Script 2_reduce_missing.py will end up creating some empty files or alignments with 
less than 3 sequences. This script will just grab the probe sequence for that locus and create 
a new fasta that can be used for the consensensus creation to map to in step 6
'''

import os
import shutil

fh_p = open('../Salamander-input-seq.fas','r')
missing = []
headers = []
seqs = []
refs = []

#make outdir
outdir = '../missing_fasta'
if not os.path.exists(outdir):
	os.mkdir(outdir)

#write ref fastas to a 2 lists
for line in fh_p:
	line=line.strip()
	if line.startswith('>'):
		headers.append(line.strip('>'))
	else:
		seqs.append(line)
		
#combine to a single list
print 'combining fasta file'
for h,s in zip(headers,seqs):
	refs.append(h+'\t'+s)

#write missing loci to a missing list
for fasta in os.listdir('.'):
	if fasta.endswith('fasta') and fasta.startswith('set'):
		fh = open(fasta,'r')
		if len(fh.readlines()) < 4:
			print len(fh.readlines())
			missing.append(fasta.split('.')[0])
		elif len(fh.readlines()) == 4:
			print len(fh.readlines())
			missing.append(fasta.split('.')[0])
			
#compare the two lists and write the seqs to new files with 3 copies as before
for seq in refs:
	for locus in missing:
		if seq.split('\t')[0].split('|')[0] == locus or seq.split('\t')[0].replace('|','-'):
			print locus
			newname = locus.replace('|','-')+'.fasta'
			fh_temp = open(newname,'a')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1]+'\n')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1]+'\n')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1])
			fh_temp.close()
			shutil.move(newname,outdir)

			