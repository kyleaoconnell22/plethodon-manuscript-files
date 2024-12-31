import os
import sys

'''
python sum_bams.py bam_dir_path reference_path

'''

reference = sys.argv[1] #reference path (root name, needs to be indexed first)
min_cov = 5

def bcf_cons(pileup,out_fasta,cov):
    with open(pileup) as f:
        content = [x.strip('\n') for x in f.readlines()]

    loci = []
    seq_dict = {}
    for line in content:
        if '#' not in line:
            # Split the tab delimited lines into their segments
            element = line.split('\t')
            # Create a list of all loci names
            seq_name = element[0]
            if seq_name not in loci:
                loci.append(seq_name)
            # By default call every position a uncertainty
            basecall = "N"

            # Turn all lower case values in upper case
            sample = element[4].upper()
            # make a directory with all different basecalls and count their occurences
            calls = dict((letter,sample.count(letter)) for letter in set(sample))

            # The basecall in the reference
            reference = element[2]
            # These characters signal a match with the reference
            match_ref = "." ","
            # List of base characters
            bases = "A" "G" "C" "T"

            # find out how many agreements with reference are among the basecalls. These are all . and , basecalls listed in match_ref.
            # reset the counter before every round (every line in file aka position in sequence)
            list_matches = 0
            for key,value in calls.items():
                if key in match_ref:
                    list_matches += value
            if list_matches >= cov:
                basecall = reference

            # find if there are any well supported SNPs and make the most prominent call
            for key in sorted(calls, key=calls.get, reverse=True):
                if key in bases:
                    if int(calls[key]) >= cov:
                        if int(calls[key]) >= list_matches:
                            basecall = key
                            break
            # add the final basecall to the dictionary and to the respective key if it already exists, otherwise create new key
            seq_dict.setdefault(element[0],[])
            seq_dict[element[0]].append(basecall)
    # Join all basecalls for each key (=locus) into one sequence and deposit in new dictionary
    concat_basecalls = {}
    for key, value in seq_dict.items():
        concat_basecalls[key] = "".join(value)
    #print concat_basecalls

    with open(out_fasta, "w") as f:
        for k, v in concat_basecalls.items():
            f.write(">" + k+ "\n")
            f.write(v+ "\n")
    return out_fasta


for bam_file in os.listdir('.'):
	if bam_file.endswith('.bam'):
		name_base=bam_file.split('.')[0]
		mp_out = bam_file.split('.')[0]+'.bcf'
		#print 'running mpileup on ', bam_file
		#mp_cmd = 'bcftools mpileup -A -f %s %s > %s'%(reference,bam_file,mp_out)
		#os.system(mp_cmd) #run it
		
		print 'finding consensus for ', name_base
		fasta_file = name_base + '_temp.fasta'
		fasta_file = bcf_cons(mp_out,fasta_file,min_cov)


'''
		mpileup_file = name_base + '.bcf'
		out = open(mpileup_file, 'w') #figure this out!!!!
		mp = subprocess.Popen(mpileup, stdout=out, stderr=subprocess.PIPE)
		mp.communicate()
		mp.wait()
		
		fasta_file = name_base  + '_temp.fasta'
		fasta_cons = bcf_cons(out,fasta_file,min_cov)

#>set6-AMEXTC-0340000064262-ndnf-exon-3_cinereus_IGG_Blood_64114 |set6-AMEXTC-0340000064262-ndnf-exon-3
out_fasta = open('joined_unphased_fastas.fasta','w')
for file in os.listdir(path):
	if file.endswith('_temp.fasta'):
		print file
		s = file.split('.')[1]
		fh = open(file,'r')
		for line in fh:
			if line.startswith('>'):
				out_fasta.write(line+'_'+s + " |" + line.strip('>'))
			else:
				out_fasta.write(line)

'''