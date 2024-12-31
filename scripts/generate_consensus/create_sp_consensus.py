import sys
import Bio
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
import os
import subprocess as sp
import shutil
import argparse
import time
'''
#this is to create species specific consensus for mapping probes to a a reference
Make sure you run : source activate secapr_env
python create_sp_consensus.py in_dir species_list.txt
species list is just list of target species to iterate through
Will sort alignment output from secapr aln dir by species, find consensus for that species, then write
that to out file with name species_consensus.fasta
If all loci are not present, you can use the next script to fill in missing sequences: fill_gaps_consensus.py
Kyle O'Connell
June 2020
'''
in_dir = sys.argv[1] 
fh_sp = open(sys.argv[2],'r') #needs to be species list
sp = []
os.chdir(in_dir) #change into target dir
for line in fh_sp:
	line=line.strip()
	sp.append(line)
for species in sp:
	out_path = species+'_consensus.fasta'
	if not os.path.exists(out_path):
		fh_out = open(out_path,'a')


#remove empty fasta if they exist
for fasta in os.listdir('.'):
	if fasta.endswith('.fasta') and not fasta.endswith('consensus.fasta'):
		if len(open(fasta).readlines()) < 2:
			os.remove(fasta)

#testing cons module from secapr_env
#now create species specific consensus
for fasta in os.listdir('.'): #iterate through
	if fasta.endswith('.fasta') and not fasta.endswith('consensus.fasta'): #grab only fasta files
		alignment = AlignIO.read(fasta, 'fasta') #read fasta to alignIO object
		for species in sp: #take one species at a time
			temp_seq = species+'_'+fasta.split('.')[0]+'_temp.fasta'
			temp_list=[]
			fh_temp = open(temp_seq,'a')
			for sequence in alignment: #take one sequence at a time in the alignment obj
				if not 'FFS' in sequence.id: #exclude FFS for consensus
					if sequence.id.split('_')[0] == species:
						fh_temp.write('>'+sequence.id+'\n'+str(sequence.seq)+'\n')
						temp_list.append(sequence.id)

			#if empty file, skip it
			if len(open(temp_seq).readlines()) <= 1:
				pass
			else: #generate consensus
				con_file = 'con.fasta'
				header = fasta.split('.')[0]
				cons_cmd = "cons -sequence %s -outseq %s -name %s -plurality 0.1 -setcase 0.1" %(temp_seq,con_file,header)
				os.system(cons_cmd)
				con = SeqIO.read(con_file,'fasta')
				if len(con)>100:
					fh_con = open(con_file,'r')
					out_path = species+'_consensus.fasta'
					fh_out = open(out_path,'a')
					for line in fh_con:
						fh_out.write(line)
			os.remove(temp_seq)

#module for testing the similarity between species sets
#writes the sequence for each species of the target locus hard coded to a single file
target = '2310' #total number of targets
fh_out=open('sanity_check.fasta','a')
for filetype in os.listdir('.'):
	if filetype.endswith('consensus.fasta'):
		with open(filetype, "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				if record.id == target:
					fh_out.write('>'+filetype.split('_')[0]+'_')
					SeqIO.write(record,fh_out,'fasta')

#this just rewrites the file above as single line
fh_sanity = open('sanity_check.fasta','r')
fh2 = open(target+'_consensus.fasta','a')
for line in fh_sanity:
	if line.startswith('>'):
		line=line.strip()
		fh2.write('>'+line.split('_')[0].split('>')[1]+'_'+target+'\n')
	else:
		fh2.write(line)
os.remove('sanity_check.fasta')	

