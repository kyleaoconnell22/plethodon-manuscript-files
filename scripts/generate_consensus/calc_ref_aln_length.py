import os
import sys
from Bio import SeqIO
import numpy

for fas in os.listdir('.'):
	seq_len = []
	name = fas.split('.')[0]
	ref = SeqIO.parse(fas,'fasta')
	for record in ref:
		seq_len.append(len(record.seq))
		#print len(record.seq)
	print name, numpy.mean(seq_len)
		
