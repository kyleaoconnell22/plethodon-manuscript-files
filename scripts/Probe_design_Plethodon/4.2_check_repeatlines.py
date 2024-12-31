fh = open("4.3-reduced_final_seqs.fa", "r")
fh_out = open("repeat_lines", 'a')
seqlist = []
seqdict = {}
for line in fh:
	line = line.strip()
	if line.startswith(">"):
		seqdict[line]=''
		if line not in seqlist:
			seqlist.append(line)
		else:
			fh_out.write(line+'\n')