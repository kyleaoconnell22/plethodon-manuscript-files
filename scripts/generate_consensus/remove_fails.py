import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
import os
import subprocess as sp
import shutil
import argparse
import random

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-f", "--in_fasta", required=True, help="REQUIRED: input fasta aln")
    return parser.parse_args()

def remove_fails(in_dir,in_fasta):
	os.chdir(in_dir)
	
	new_out = in_fasta.split('_')[0]+'_consensus_trimmed_cleaned.fasta'
	fh_out=open(new_out,'a')
	alignment = AlignIO.read(in_fasta, 'fasta')
	
	for item in alignment:
		seq = str(item.seq)
		dash = item.seq.count('-')
		N = item.seq.count('N')
		miss=int(dash)+int(N)
		leng = len(item.seq)
		perc = float(miss)/(len(item.seq))
		#round to two decimal places
		perc = round(perc,2)
		print item.id, perc
		if perc < 0.75:
			fh_out.write(str('>'+item.id+'\n'+seq+'\n'))
		else:
			print item.id
def main():
	#define the arguments
	args = get_args()
	remove_fails(args.in_dir,args.in_fasta)
	
if __name__ == '__main__':
    main()	
