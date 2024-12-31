import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
import os
import subprocess as sp
import shutil
import argparse
import random
from Bio.Nexus import Nexus
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

'''
two modules here
first will convert fasta to nex
Second will concatenate dir of nex files and convert to phy for phylogenetic analysis
Can remove nex files if desired

'''
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-o", "--out_dir", required=True, help="REQUIRED: The full path to the dir with the nexus files")
    parser.add_argument("-p", "--prefix", required=True, help="REQUIRED: prefix for the final concatentated file")
    return parser.parse_args()

def create_out(out_dir):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
#convert all fastas to nex using SeqIO and then move the nex files to out_dir
def convert_fasta(in_dir,out_dir):
	os.chdir(in_dir)
	for file in os.listdir('.'):
		if file.endswith('.fasta'):
			input_handle = open(file, "rU")
			out=file.split('.')[0]+'.nex'
			output_handle = open(out, "w")
			SeqIO.convert(file, "fasta", out, "nexus", generic_dna)

			shutil.move(out,out_dir)

#create file list of all nex files
#use Nexus.combine to concat all nex files
#then use SeqIO to convert to phy for raxml or iqtree
#then remove all nex files, this line can be commented out if this is not desired	
def concatenate(out_dir,prefix):
	os.chdir(out_dir)
	file_list = []
	for file in os.listdir('.'):
		if file.endswith('.nex'):
			file_list.append(file)
	nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]

	combined = Nexus.combine(nexi)

	cat = prefix+'.nex'
	cat_out = open(prefix+'.phy',"w")

	combined.write_nexus_data(filename=open(cat, 'w'))

	for file in os.listdir('.'):
		if file.endswith('single.nex'):
			os.remove(file)

def convert_final(out_dir,prefix):
	os.chdir(out_dir)
	fh_in = open(prefix+'.nex','r')
	fh_out = open(prefix+'.phy','a')
	
	for line in fh_in:
		line=line.strip()
		if line.startswith('#'):
			pass
		elif line.startswith('begin') or line.startswith('format') or line.startswith('matrix'):
			pass
		elif line.startswith('dimensions'):
			ntax = line.split(' ')[1].split('=')[1]
			nchar = line.split(' ')[2].split('=')[1].strip(';')
			print ntax
			print nchar
			fh_out.write(ntax+' '+nchar+'\n')
		elif line.startswith('matrix'):
			pass
		elif line.startswith('char'):
			pass
		elif line.startswith('end'):
			pass
		elif line.startswith(';'):
			pass
		else:
			print line
			fh_out.write(line.replace('?','N')+'\n')
	
def main():
	#define the arguments
	args = get_args()
	create_out(args.out_dir)
	convert_fasta(args.in_dir,args.out_dir)
	concatenate(args.out_dir,args.prefix)
	convert_final(args.out_dir,args.prefix)

if __name__ == '__main__':
    main()	
