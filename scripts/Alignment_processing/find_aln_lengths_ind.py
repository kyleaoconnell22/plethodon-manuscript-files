import argparse
import os
import sys
import subprocess as sp
import shutil
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np
import os
import sys
fh = open(sys.argv[1],'r') #sample list
samples = []
for line in fh:
	line=line.strip()
	samples.append(line)
fh_out=open('../ind_lens_summary.txt','a')
fh_out.write('sample'+'\t'+'ave_len'+'\n')
for sample in samples:
	l = [] 
	for filetype in os.listdir('.'):
		if filetype.endswith('.fasta'):
			with open(filetype, "rU") as handle:
				for record in SeqIO.parse(handle, "fasta"):
					if record.id == sample:
						lseq = len(record.seq.strip('-'))
						if lseq>0:
							l.append(len(record.seq.strip('-')))
	avel = round(float(np.mean(l)),2)
	print str(avel), sample
	fh_out.write({0}+'\t'+{1}+'\n'.format(sample,str(avel)))
	