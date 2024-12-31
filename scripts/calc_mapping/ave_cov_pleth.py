'''
In R, run the following to get the input for this script
cov = read.table('average_cov_per_locus.txt',header=T)
covT=t(cov)
write.csv(covT,'covT.csv')
Then you need to get rid of the " " manually in bbedit

'''

fr = []
bl = []
ffs = []

import numpy
fh = open('covT.csv','r')
fh_out = open('ave_cov.txt','a')
fh_out.write('Sample'+'\t'+'ave_cov'+'\t'+'data_type'+'\n')
for line in fh:
	if line.startswith('Samples'):
		pass
	elif line.startswith('locus'):
		pass
	else:
		tmp_cov = []
		line=line.strip()
		s = line.split(',')[0]
		c = line.split(',')[1:]
		for item in c:
			if round(float(item),2) > 0.1:
				tmp_cov.append(float(item))
		ave_cov = round(float(numpy.mean(tmp_cov)),2)
		fh_out.write(s+'\t'+str(ave_cov)+'\t'+s.split('_')[2]+'\n')
		if 'Fresh' in s:
			fr.append(ave_cov)
		elif 'Blood' in s:
			bl.append(ave_cov)
		elif 'FFS' in s:
			ffs.append(ave_cov)

ave_fr = str(round(numpy.average(fr),2))
ave_bl = str(round(numpy.average(bl),2))
ave_ffs = str(round(numpy.average(ffs),2))

sd_fr = str(round(numpy.std(fr),2))
sd_bl = str(round(numpy.std(bl),2))
sd_ffs = str(round(numpy.std(ffs),2))

min_fr = str(round(min(fr),2))
max_fr = str(round(max(fr),2))
min_bl = str(round(min(bl),2))
max_bl = str(round(max(bl),2))
min_ffs = str(round(min(ffs),2))
max_ffs = str(round(max(ffs),2))



fh_out.write('\n')
fh_out.write('ave cov fresh'+'\t'+ave_fr+ '\t'+ 'range' + '\t' + min_fr +'-'+ max_fr + '\t' + 'SD = ' + sd_fr + '\n')
fh_out.write('avecov blood'+'\t'+ave_bl+ '\t'+ 'range' + '\t' + min_bl +'-'+max_bl +  '\t' + 'SD = ' + sd_bl + '\n')
fh_out.write('ave cov ffs'+'\t'+ave_ffs+'\t'+ 'range' + '\t'+ min_ffs +'-'+max_ffs + '\t' + 'SD = ' + sd_ffs +'\n')
