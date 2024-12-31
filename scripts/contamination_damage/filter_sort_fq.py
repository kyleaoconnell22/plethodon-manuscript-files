from Bio import SeqIO
import os

'''
script is designed to fix the mess caused by removing exogenous contamination reads using 
fastq_screen, which analyzes the R1 and R2 separately. So it will find the reads not shared by the two 
files and then filter them out and write new filtered files to a dir called sample_filt
that have exactly the same number of reads
This is needed for mapping with bwa or bowtie2

Kyle O'Connell
Oct. 2020
kyleaoconnell22@gmail.com

'''

#iterate through list dir
for dir in os.listdir('.'):
	if os.path.isdir(dir):  #make sure its dir not file
		os.chdir(dir) #change into each sample dir (could skip this and modify paths below)
		print '\n'
		print 'working with', dir
		R1_rec = [] #lists for counting reads in each file
		R2_rec = []

		#chunk for finding differences between the two files
		for handle in os.listdir('.'):
			if handle.endswith('READ1.fastq'): #R1 file (ignores single reads files if present)
				for record in SeqIO.parse(handle, "fastq"):
					R1_rec.append(record.id) #just write the read headers to a list
			elif handle.endswith('READ2.fastq'):
				for record in SeqIO.parse(handle, "fastq"):
					R2_rec.append(record.id)
		print 'len R1 file', len(R1_rec) #print number of reads in R1 and R2 (diff if filtered by fastq_screen)
		print 'len R2 file', len(R2_rec)
		
		#find reads not shared by the two files
		m1 = list((set(R1_rec)-set(R2_rec)))
		m2 = list((set(R2_rec)-set(R1_rec)))
		missing = set(m1+m2)
		print 'diff between the R1 and R2 = ', len(missing)	
	
		#chunk for filtering missing ones and writing to new files
		out_dir = '../'+dir+'_filt' #file name + filt
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)	
		
		outR1 = out_dir+'/'+ dir.split('.')[0]+'_clean-READ1.fastq' #out files
		outR2 = out_dir+'/'+ dir.split('.')[0]+'_clean-READ2.fastq'
		R1_keep = [] #shared reads
		R2_keep = []

		#write reads to keep lists if shared between R1 and R2
		for handle in os.listdir('.'):
			if handle.endswith('READ1.fastq'):
				for record in SeqIO.parse(handle, "fastq"):
					if record.id not in missing:
						R1_keep.append(record) 
			elif handle.endswith('READ2.fastq'):
				for record in SeqIO.parse(handle, "fastq"):
					if record.id not in missing:
						R2_keep.append(record)
	
		#write keep lists to outfiles in out dir
		print 'len R1_keep', len(R1_keep)
		print 'len R2_keep', len(R2_keep)
		SeqIO.write(R1_keep, outR1, "fastq")
		SeqIO.write(R2_keep, outR2, "fastq")
		os.chdir('../') #make sure to change back out to base dir