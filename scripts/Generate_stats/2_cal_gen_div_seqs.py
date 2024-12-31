import dendropy
from dendropy.calculate import popgenstat
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

'''
need path to species list, which is just hte species you want to iterate through to calc stats
Assumes you have ffs, blood, and fresh

Assumes you are in the listdir already

'''

specieslist = []

fh_sp=open(sys.argv[1],'r')
for line in fh_sp:
	line=line.strip()
	specieslist.append(line)

for species in specieslist:
	ffs =[]
	bl = []
	fr = []
	
	for file in os.listdir('.'):
		if file.endswith('.fasta'):
			aln = AlignIO.read(file, 'fasta')
			temp_ffs = open('ffs.temp.fasta','a')
			temp_bl = open('bl.temp.fasta','a')
			temp_fr = open('fr.temp.fasta','a')
			seqs = []
			for record in aln:
				if record.id not in seqs:
					seqs.append(record.id)
					if record.id.startswith(species):
						if 'FFS' in record.id:
							SeqIO.write(record,temp_ffs,'fasta')
						elif 'Blood' in record.id:
							SeqIO.write(record,temp_bl,'fasta')
						elif 'Fresh' in record.id:
							SeqIO.write(record,temp_fr,'fasta')

			alnffs = dendropy.DnaCharacterMatrix.get(file=open("ffs.temp.fasta"), schema="fasta")
			os.remove('ffs.temp.fasta')
			os.remove('bl.temp.fasta')
			os.remove('fr.temp.fasta')

			for seq in alnffs.values():
				print file, len(seq)
				
			####********something is going on with the last sequence, need to figure out why that is! ####
'''
    		alnbl = dendropy.DnaCharacterMatrix.get(path="bl.temp.fasta", schema="fasta")
			alnfr = dendropy.DnaCharacterMatrix.get(path="fr.temp.fasta", schema="fasta")
			if len(alnffs)> 1:
				ffs_pi = dendropy.calculate.popgenstat.nucleotide_diversity(alnffs, ignore_uncertain=False)
			
			
			
						
			alnffs = seqs = dendropy.DnaCharacterMatrix.get(path="ffs.temp.fasta", schema="fasta")
			alnbl = seqs = dendropy.DnaCharacterMatrix.get(path="bl.temp.fasta", schema="fasta")
			alnfr = seqs = dendropy.DnaCharacterMatrix.get(path="fr.temp.fasta", schema="fasta")
			print alnffs

			heFFS = dendopy he of alnffs
			bl and fr same
			
			ffs.append(heFFS)
			bl.append(hebl)
			fr.append(hefr)
	
	#calc averages of he across each of the three lists
	if len(ffs)>0:
		ave_he_ffs = np.mean(ffs)
		fh_out.write(Heterozygosity of {0} FFS = {1} +'\n'.format(species,ave_he_ffs))
	elif len(bl)>0:
		ave = xxx
		fh_out.write(Heterozygosity of {0} Blood = {1} +'\n'.format(species,ave_he_bl))
	elif len(fr)>0:
		ave = XX
		fh_out.write(Heterozygosity of {0} Blood = {1} +'\n'.format(species,ave_he_bl))
'''