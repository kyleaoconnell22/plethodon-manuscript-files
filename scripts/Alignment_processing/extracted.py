fh=open('extracted_target_contigs_all_samples.fasta','r')
samples = []
for line in fh:
	line=line.strip()
	if line.startswith('>'):
		sample=line.split('>')[1].split(' |')[0].split('_')[-1]
		if sample not in samples:
			samples.append(sample)
			print sample

print len(samples)
for item in samples:
	print item