import os
import sys
import numpy
fh_out = open('mapping_stats_summary.txt','w')
fh_out.write('sample'+'\t'+'mapping_percent'+'\n')
fr = []
bl = []
ffs = []

for filetype in os.listdir('.'):
	if filetype.endswith('.stats'):
		s = filetype.split('.')[0]
		fh=open(filetype,'r')
		line = fh.readlines()[2]
		map = float(line.split('(')[1].split('-')[0].strip('%:'))
		if 'Blood' in filetype:
			bl.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')
		elif 'Fresh' in filetype:
			fr.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')
		elif 'FFS' in filetype:
			ffs.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')

ave_fr = str(round(numpy.average(fr),2))
ave_bl = str(round(numpy.average(bl),2))
ave_ffs = str(round(numpy.average(ffs),2))

min_fr = str(round(min(fr),2))
max_fr = str(round(max(fr),2))
min_bl = str(round(min(bl),2))
max_bl = str(round(max(bl),2))
min_ffs = str(round(min(ffs),2))
max_ffs = str(round(max(ffs),2))



fh_out.write('\n')
fh_out.write('ave_mapping_fr'+'\t'+ave_fr+ '\t'+ 'range' + '\t' + min_fr +'-'+ max_fr + '\n')
fh_out.write('ave_mapping_bl'+'\t'+ave_bl+ '\t'+ 'range' + '\t' + min_bl +'-'+max_bl + '\n')
fh_out.write('ave_mapping_ffs'+'\t'+ave_ffs+'\t'+ 'range' + '\t'+ min_ffs +'-'+max_ffs +'\n')

