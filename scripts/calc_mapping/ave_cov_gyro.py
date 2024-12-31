import numpy as np
'''
In R, run the following to get the input for this script
The input for the R script comes from secapr after the mapping step (s6)
cov = read.table('average_cov_per_locus.txt',header=T)
covT=t(cov)
write.csv(covT,'covT.csv')
Then you need to get rid of the " " manually in bbedit

'''
ffs = []
al = []
fr = []
rad = []

fh = open('covT.csv','r')
fh_out = open('ave_cov.txt','w')
fh_out.write('Sample'+'\t'+'ave_cov'+'\t'+'sd'+'\t'+'data_type'+'\n')
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
				if s.split('_')[1] == 'FFS':
					ffs.append(float(item))
				elif s.split('_')[1] == 'allozyme':
					al.append(float(item))
				elif s.split('_')[1] == 'fresh':
					fr.append(float(item))
				elif s.split('_')[1] == 'rad':
					rad.append(float(item))
		#ind level stats	
		ave_cov = str(round(np.mean(tmp_cov),2))
		sd_cov = str(round(np.std(tmp_cov),2))
		#print s, ave_cov
		fh_out.write(s+'\t'+str(ave_cov)+'\t'+ '\t'+ sd_cov +'\t'+ s.split('_')[1]+'\n')
		
#data type level stats
fh_out.write('\n'*2)
ave_ffs = str(round(np.mean(ffs),2))
sd_ffs = str(round(np.std(ffs),2))
ave_al = str(round(np.mean(al),2))
sd_al = str(round(np.std(al),2))
ave_fr = str(round(np.mean(fr),2))
sd_fr = str(round(np.std(fr),2))
ave_rad = str(round(np.mean(rad),2))
sd_rad = str(round(np.std(rad),2))

#write to out file
fh_out.write('ave_ffs = {0}, sd = {1}'.format(ave_ffs,sd_ffs)+'\n')
fh_out.write('ave_al = {0}, sd = {1}'.format(ave_al,sd_al)+'\n')
fh_out.write('ave_fr = {0}, sd = {1}'.format(ave_fr,sd_fr)+'\n')
fh_out.write('ave_rad = {0}, sd = {1}'.format(ave_rad,sd_rad)+'\n')

