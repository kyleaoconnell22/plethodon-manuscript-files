import os
import numpy

'''
first need to cat R1 and R2 contam, so that the screen file has values for both
then I write these values to a list and find the mean
Write this to the outfile

If only use r1 then need to get rid of the mean part and just write the first value to fh_out

#samples is like:
sample1_ffs
sample1_allozyme
sample1_fresh
sample1_rad


for x in *_clean-READ1_screen.txt; do cat $x_clean-READ1_screen.txt $x_clean-READ2_screen.txt > $x_comb_screen.txt; done 
This wasnt working so I used regular expression in text wrangler to manually cat them
'''

fh_out = open('out_contam.txt','a')
fh_out.write('sample'+'\t'+ 'total_reads'+'\t'+'%total_contam'+'\n')

all_values = []
FFS = []
Blood = []
Fresh = []
RAD = []

for filetype in os.listdir('.'):
	if '_comb_screen.txt' in filetype: #only get combined screen out files (includes R1 and R2)
		cont = [] #empty list for R1 and R2 value
		reads_t = []
		s = filetype.split('_')[0]+'_'+filetype.split('_')[1]+'_'+filetype.split('_')[2]\
			+'_'+filetype.split('_')[3] #sample name
		fh = open(filetype,'r')
		for line in fh:
			line=line.strip()
			if line.startswith('Human'):
				reads = line.split('\t')[1]
				reads_t.append(int(reads))
			if line.startswith('%'): #just grab %unmapped
				unmapped = float(line.split(':')[1].strip(' ')) #value not mapped to human
				contam = float(100-unmapped) #find value of human mapping
				cont.append(contam) #write R1 and R2 to list
		#print cont
		both_reads = str(round(numpy.mean(cont),2)) #find mean
		#print both_reads
		sum_reads = sum(reads_t)
		fh_out.write(s + '\t' + str(sum_reads) + '\t' + both_reads + '\n') #write out
		out = s + ',' + both_reads
		all_values.append(out)
		

for i in all_values:
	if 'FFS' in i:
		FFS.append(float(i.split(',')[1]))
	elif 'Blood' in i or 'allozyme' in i:
		Blood.append(float(i.split(',')[1]))
	elif 'Fresh' in i or 'fresh' in i:
		Fresh.append(float(i.split(',')[1]))
	elif 'rad' in i:
		RAD.append(float(i.split(',')[1]))

ave_FFS = str(round(numpy.mean(FFS),2)) #find mean
min_FFS = str(round(min(FFS),2)) #find min
max_FFS = str(round(max(FFS),2)) #find max

ave_Blood = str(round(numpy.mean(Blood),2)) #find mean
min_Blood = str(round(min(Blood),2))
max_Blood = str(round(max(Blood),2))


ave_Fresh = str(round(numpy.mean(Fresh),2)) #find mean
min_Fresh = str(round(min(Fresh),2))
max_Fresh = str(round(max(Fresh),2))

if len(RAD) > 2:
	ave_RAD = str(round(numpy.mean(RAD),2)) #find mean
	min_RAD = str(round(min(RAD),2))
	max_RAD = str(round(max(RAD),2))


#write summary values to fh_out
fh_out.write('\n' + 'average contamination FFS = ' + ave_FFS + '\n')
fh_out.write('range of FFS = ' + min_FFS + '-' + max_FFS + '\n')

fh_out.write('average contamination Blood/Allozyme = ' + ave_Blood + '\n')
fh_out.write('range of Blood/Allozyme = ' + min_Blood + '-' + max_Blood + '\n')

fh_out.write('average contamination Fresh = ' + ave_Fresh + '\n')
fh_out.write('range of Fresh = ' + min_Fresh + '-' + max_Fresh + '\n')

if len(RAD) > 2:
	fh_out.write('average contamination RAD = ' + ave_RAD + '\n')
	fh_out.write('range of RAD = ' + min_RAD + '-' + max_RAD + '\n')

	
		
		