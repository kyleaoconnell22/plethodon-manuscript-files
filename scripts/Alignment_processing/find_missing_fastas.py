'''
script will create a fasta file from the input probe locus file to include in the locus mapping
Needs a list of loci that secapr excludes because they have less than 3 samples
'''

fh = open('missing_list.txt','r')
infile = '4_target_contigs/pleth_bestsamples/reference_fasta_header_info.txt'
fh2 = open(infile,'r')
fh_p = open('fasta_ref/Salamander-input-seq.fas','r')

missing = []
all_loci = []
missing_new = []
headers = []
seqs = []
x = []

for line in fh:
	line=line.strip()
	missing.append(line)
for line in fh2:
	line=line.strip()
	all_loci.append(line)

for line in fh_p:
	line=line.strip()
	if line.startswith('>'):
		headers.append(line.strip('>'))
	else:
		seqs.append(line)

print 'renaming missing loci'		
for l in all_loci:
	for locus in missing:
		if l.split('\t')[0] == locus:
			missing_new.append(l.split('\t')[1])

print 'combining fasta file'
for h,s in zip(headers,seqs):
	x.append(h+'\t'+s)

for seq in x:
	for locus in missing_new:
		if seq.split('\t')[0] == locus:
			newname = locus.replace('|','-')+'.fasta'
			fh_temp = open(newname,'a')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1]+'\n')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1]+'\n')
			fh_temp.write('>'+seq.split('\t')[0]+'\n'+seq.split('\t')[1])

			
			