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

#assign base dir and change to it, I recommend Trimmed
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#Reference Genome spliced into Exons
filename = "GCF_000935625.1_ASM93562v1_genomic.fna"
fh = open(filename, 'r')

#transcript file
filename2 = "/home/oconnellk/Desktop/exon_capture/probe_fastas/RHACO-9715-seq_kebi.fasta"
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
proc_blastn = sp.Popen(['blastn', '-db', db, '-query', query, '-out', outfile, '-num_threads', '4', '-outfmt',
'6 sseqid qseqid evalue bitscore qstart qend'])
proc_blastn.wait()

#stringent, does not allow multiple hits (typically)
#proc_blastn = sp.Popen(['blastn', '-db', db, '-query', query, '-out', outfile, '-num_threads', '4', '-outfmt',
#'6 sseqid qseqid evalue bitscore qstart qend', '-evalue', '1e-10','-best_hit_score_edge', '0.05','-best_hit_overhang', '0.25'])
#proc_blastn.wait()


fh.close()
fh2.close()
