''''
Usage: python 3_blast_Transcripts_to_Ref [PseudoRefGenome.fa}] [Transcript_CDS_file.fa]

1) filters CSV file down to one exon per locus

------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
Sept 2018
------------------------
'''

#________________________________________#
import sys
import os
import subprocess as sp
#________________________________________#
print '\n'
print "###############################################"
print "###############################################"

#assign base dir and change to it, I recommend Trimmed
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#Pseudo reference genomes (only exons)
filename = "vision_seqs.fasta"
fh = open(filename, 'r')

#transcript file
filename2 = "2.1_spliced_reference_exons.fa"
fh2 = open(filename2, 'r')

##########################################################
#create Blast DB of ref genome
#subject
proc_blastdb = sp.Popen(['makeblastdb', '-in', filename, '-dbtype', 'nucl'])
proc_blastdb.wait()
db = filename

#query is transcriptome sequences file
query = filename2	
outfile = query + '.blastout'

#loose params (allows multiple hits)
#proc_blastn = sp.Popen(['blastn', '-db', db, '-query', query, '-out', outfile, '-evalue', '0.001', '-num_threads', '4', '-outfmt',
#'6 sseqid qseqid evalue bitscore qstart qend'])
#proc_blastn.wait()

#stringent, does not allow multiple hits (typically)
proc_blastn = sp.Popen(['tblastx', '-db', db, '-query', query, '-out', outfile, '-num_threads', '4', '-outfmt',
'6 sseqid qseqid evalue bitscore qstart qend', '-evalue', '20','-best_hit_score_edge', '0.05','-best_hit_overhang', '0.25'])
proc_blastn.wait()
i = 0
fh_blast = open(outfile, 'r')
for line in fh_blast:
	i = i + 1

print "There were {0} blast hits".format(i)
fh.close()
fh2.close()
print "###############################################"
print "###############################################"
print '\n'