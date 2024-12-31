fh = open("4_spliced_transcripts.fa", "r")
fh_out = open("4_filtered_final_sequences.fa", 'a')
seqlist = []
seqdict = {}
for line in fh:
	line = line.strip()
	if line.startswith(">"):
		name = line
		seqdict[name]=''
	else:
		seq = line
		seqdict[name]=seq
for key,value in seqdict.iteritems():
	fh_out.write(str(key)+ '\n' + str(value)+'\n')
