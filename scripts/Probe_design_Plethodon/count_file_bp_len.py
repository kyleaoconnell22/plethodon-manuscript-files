#________________________________________#
import sys
import os
import subprocess as sp
import random
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO

#________________________________________#

#final target exon fasta file
filename = sys.argv[1]
#count total number of loci
i = 0
len_list = []
for record in SeqIO.parse(filename, "fasta"):
    i = i + len(record.seq)
    bp = len(record.seq)
    len_list.append(bp)

ave_len = sum(len_list)/len(len_list)
max_len = max(len_list)
min_len = min(len_list)	
		

print "total bp included in pared down target list = ", i
print "ave length of seq in file is " , ave_len
print "total length of functional genes = 53514"
print "total allowed for phylo genes = 746486"
print "max = " , max_len
print "min = " , min_len