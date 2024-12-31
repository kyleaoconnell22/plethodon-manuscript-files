import os
import numpy

'''
first need to cat R1 and R2 contam, so that the screen file has values for both
then I write these values to a list and find the mean
Write this to the outfile

If only use r1 then need to get rid of the mean part and just write the first value to fh_out
'''

fh_out = open('out_contam.txt','a')
fh_out.write('sample'+'\t'+'%human'+'\n')
for filetype in os.listdir('.'):
	if 'screen' in filetype: #only get screen out files (includes R1 and R2)
		cont = [] #empty list for R1 and R2 value
		s = filetype.split('_')[0]+'_'+filetype.split('_')[1] #sample name
		fh = open(filetype,'r')
		for line in fh:
			line=line.strip()
			if line.startswith('Human'): #ecoli was zero, could add other values here too
				unmapped = float(line.split('\t')[3]) #value not mapped to human
				human = float(100-unmapped) #find value of human mapping
				cont.append(human) #write R1 and R2 to list
		#print cont
		both_reads = str(round(numpy.mean(cont),2)) #find mean
		#print both_reads
		fh_out.write(s + '\t' + both_reads + '\n') #write out

		