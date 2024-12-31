import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import sys

remove=[]
sample_list = []
fh_pa = open(sys.argv[1],'r')

for line in fh_pa:
	line=line.strip()
	remove.append(line)

for file in os.listdir('.'):
	if file.endswith('.fasta'):
		for record in SeqIO.parse(file, "fasta"):
			sample_list.append(record.id)

sample_set = set(sample_list)
N = len(sample_set)

for file in os.listdir('.'):
	if file.endswith('.fasta'):
		aln = AlignIO.read(file, 'fasta')
		for record in aln:
			if len(record.seq)<50:
				remove.append(file)
for file in os.listdir('.'):
	if file.endswith('.fasta'):
		#open fasta to count seqs inside
		aln = AlignIO.read(file, 'fasta')
		x = float((len(aln)/2))/float((N/2))
		if x > 0.9: #require 10% of samples to be present in every aln to keep it
			remove.append(file)
print 'seqs to remove that captured poorly: '

for aln in set(remove):
	#print aln
	os.remove(aln)

for file in os.listdir('.'):
	if file.startswith('ocus'):
		newname='l'+file.split('_gblocks')[0]+'.fasta'
		os.rename(file,newname)
	elif file.startswith('set'):
		newname=file.split('_gblocks')[0]+'.fasta'
		os.rename(file,newname)

		