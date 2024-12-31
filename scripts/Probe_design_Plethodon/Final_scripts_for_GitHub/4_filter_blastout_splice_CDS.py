'''
python 4_filter_blastout_splice_CDS [Base_Dir] Notice that you need to hard code a bunch of file names

means parameters that can be coded into the script

Edit everything within the parameters section below
------------------------
written for Python 2.7
Blake O'connell and Kyle O'Connell
Sept. 2018
------------------------
'''

#________________________________________#
import sys
import os
import subprocess as sp
import random
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import numpy
import numpy as np
#________________________________________#

#Parameters (edit these according to your system)
########################
min_length = 199
trim_length = 400
gc_low = 30
gc_high = 70
'''
target_loci is because you are limited in the number of loci (really number of bp) that can be included in the probe design
so if you have more than the target loci that survive all filtering, this script (4.3 section) will pick exons randomly from all acceptable loci
'''
target_loci = 2202 
blastoutput = "Sals_annotated.fasta.blastout"
transcriptCDS = "Sals_annotated.fasta"
final_seqs = "4_filtered_final_sequences.fa"
random_seqs = '4.3-random_final_seqs.fa'

########################

fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#this is the blastoutput file
filename = blastoutput
fh = open(filename, 'r')

# open the output file1 which is the trimmed blast output
file_out = '4_filtered_blastout'
fh_out = open(file_out, 'a')

# open the output file 2
file_out2 = '4_spliced_transcripts.fa'
fh_out2 = open(file_out2, 'a')

arrBlast = []
dicRefID = {}
dicTransID = {}

# loops file and creates the following:
#	- arrBlast: [refID, transID, Boolean Flag for Good Data, Line split at tabs into array from file]
#	- dicRefID: key - refID, value - array of all transIDs for that refID
#	- dicTransID: key - transID, value - array of all refIDs for that transID
for line in fh:
	line = line.strip('')
	line = line.split('\t')
	
	end = int(line[5])
	start = int(line[4])
	
	####### Edit this file name #######
	#skip splice coords shorter than 200 bp
	if (end - start) <= min_length:
		continue
	
	refID = line[0]
	transID = line[1]

	arrBlast.append([refID, transID, True, line])
	
	arrTemp = [transID]
	if refID in dicRefID:
		for val in dicRefID.get(refID):
			arrTemp.append(val)
	dicRefID[refID] = arrTemp
	
	arrTemp = [refID]
	if transID in dicTransID:
		for val in dicTransID.get(transID):
			arrTemp.append(val)
	dicTransID[transID] = arrTemp

fh.close()

# loop through blast array and mark data as Bad if:
#	- refID = another refID, but transIDs don't match
#	- transID = another transID, but refIDs don't match
for i in range(0, len(arrBlast)):
	refID = arrBlast[i][0]
	transID = arrBlast[i][1]
	
	if refID in dicRefID:
		for val in dicRefID.get(refID):
			if val != transID:
				arrBlast[i][2] = False
				
	if transID in dicTransID:
		for val in dicTransID.get(transID):
			if val != refID:
				arrBlast[i][2] = False
		
		
dicTemp = {}
dicBlast = {}

# loop through blast array again
# if Bad Data, skip
# if Unique transID, add corresponding line from blast file to blast dictionary
# else, add corresponding line from blast file to temp dictionary
#	temp dictionary: key - transID, value - array of lines from blast file
for i in range(0, len(arrBlast)):

	if arrBlast[i][2] == False:
		continue
		
	refID = arrBlast[i][0]
	transID = arrBlast[i][1]
	
	arrBlastLine = arrBlast[i][3]
	
	if len(dicTransID.get(transID)) > 1:
		arrTemp = [arrBlastLine]
		if transID in dicTemp:
			for val in dicTemp.get(transID):
				arrTemp.append(val)
		dicTemp[transID] = arrTemp
	else:
		dicBlast[transID] = arrBlastLine
		
# loop through temp dictionary and pick a random line from blast file for each repeated transID
for key, value in dicTemp.items():

	randVal = value[random.randint(0, len(value)-1)]
	dicBlast[key] = randVal


#write to filtered blast output
for k, v in dicBlast.items():
	chrom = str(v[0])
	trans = str(v[1])
	s = str(v[4])
	e = str(v[5])
	
	fh_out.write(chrom + ' ' + trans + ' ' + s + ' ' + e)
print '\n'
print "####################################"
# dicBlast has all unique, random, good data
###########################

filename = transcriptCDS
fh = open(filename, 'r')

# Loops through the CDS File and if the first value is in your Blast Dictionary, then it will print the values from the Blast Dictionary
# Boolean to indicate if next line is sequence to use
seq_line = False
i = 0

# Iterate through CDS File
for line in fh:
	line = line.strip('')
	line = line.split()
	
	if seq_line:
		# Reset Boolean
		seq_line = False
		
		seq = line[0]
		oldlen = len(seq)
		
		splice = seq[start:end]
	
		
		#filter for GC make sure to set the values you want
		gc_count = GC(splice)
		
		#filter out low or high GC content
		if gc_count > gc_low and gc_count < gc_high:
		
			#also trim down to 400, can be set to any value
			if len(splice) > trim_length:
				xlen = trim_length
				fh_out2.write('>' + str(RefID) + '\t' + str(TransID) +'\t' + str(oldlen) + '\t' + str(xlen) + '\n' + splice[0:400] + '\n')
				i = i + 1
			else:
				fh_out2.write('>' + str(RefID) + '\t' + str(TransID) +'\t' + str(oldlen) + '\t' + str(xlen) + '\n' + splice + '\n')
				i = i + 1

	# remove '>' char and check if key is in dictionary
	if line[0][1:] in dicBlast:
		# ID in Dictionary, so set values for next iteration
		values = dicBlast.get(line[0][1:])
		RefID = values[0]
		TransID = values[1]
		Eval = values[2]
		start = int(values[4])-1
		end = int(values[5])-1
		xlen = end - start
		
		# Set Boolean
		seq_line = True
		
print "{0} sequences passed gc and length filters".format(i)

# Close CDS file
fh.close()

####################################### 
#script 4.1
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
j = 0
for key,value in seqdict.iteritems():
	fh_out.write(str(key)+ '\n' + str(value)+'\n')
	j = j + 1
	
print "{0} sequences passed repeat filters".format(j)
print "\n"
####################################### 
#script 4.3 Randomly select desired number of loci

filename = final_seqs
fh = open(filename,'r')
outfile = random_seqs
fh_out = open(outfile,'a')

#intiate empy dictionary, we'll store the number of the marker as the key and all info from the list above as the value
exon_dict = {}
exon_list = []
filt_list = []
filt_dict = {}

#read in script 4 output of 2K + genes and write the seqs to a dict (exon_dict)
#I used this to capture the whole header, rather than SeqIO which just gets the seq name next to the >
for line in fh:
	line = line.strip()
	if line.startswith(">"):
		header = line
		exon_dict[header]=''
	else:
		seq = line
		exon_dict[header]=seq

#iterate through the dict and write the header and seq to a single line sep by "," and append to a list for the random sampling part		
for key,value in exon_dict.iteritems():
	if "LOC" not in key and "hypothetical" not in key and "putative" not in key:
		comb = key + "," + value
		exon_list.append(comb)


rand_list = random.sample(exon_list, target_loci)

#now write the parts of the final list to file
for item in rand_list:
	#convert to string to use .split
	line = str(item)
	#gene name
	gene = line.split(",")[0]
	#seq 
	seq = line.split(",")[1]
	#replace tabs with spaces
	gene = gene.replace("\t"," ")
	#write to output file
	fh_out.write(gene + "\n" + seq + '\n')

#count total number of loci
i = 0
seqs = []
for record in SeqIO.parse(outfile, "fasta"):
	i = i + len(record.seq)
	seqs.append(len(record.seq))
ave_len = sum(seqs)/len(seqs)
	
#Final length check
print "total bp included in pared down target list = ", i
print "ave length of seq in file is " , ave_len
print "total length of functional genes = 53514"
print "total allowed for phylo genes = 746486"
x = 0
subtr = 0
overage = i - 746486
if overage > 0:
	x = x + overage
	subtr = overage/339
else:
	x = x + 0

print y
print "you {0} bp over the limit and should remove {1} seqs from your list".format(overage,subtr)
	
fh.close()
fh_out.close()

#######################################
#script 4.2 final repeat check
fh = open("4.3-reduced_final_seqs.fa", "r")
fh_out = open("repeat_lines", 'a')
seqlist = []
seqdict = {}
rep = 0
for line in fh:
	line = line.strip()
	if line.startswith(">"):
		seqdict[line]=''
		if line not in seqlist:
			seqlist.append(line)
		else:
			fh_out.write(line+'\n')
			rep = rep + 1
print "the final file has {0} repeat lines".format(rep)
print "####################################"
print '\n'


