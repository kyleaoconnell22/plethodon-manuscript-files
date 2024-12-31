import os
'''
rename files for trimmomatic
I need to take the file names from 
glutinosus-SG-Fresh-RCB1049_R2.fastq
525132-allozyme_R2.fastq                
525132-fresh_R1.fastq                   
525132-fresh_R2.fastq                   
montanus-IGG-Fresh-KAO209_R1.fastq

Change to MoIggFrKAO209_R1.fasta

'''

for file in os.listdir('.'):
	new=file.replace('-','_')
	os.rename(file,new)
	