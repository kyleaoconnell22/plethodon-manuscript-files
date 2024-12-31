import os
import sys
import numpy
fh_out = open('mapping_stats_summary.txt','w')
fh_out.write('sample'+'\t'+'mapping_percent'+ '\t'+'SD' + '\n')
fr = []
al = []
ffs = []
rad = []

for filetype in os.listdir('.'):
	if filetype.endswith('.stats'):
		s = filetype.split('.')[0]
		fh=open(filetype,'r')
		line = fh.readlines()[2]
		map = float(line.split('(')[1].split('-')[0].strip('%:'))
		if 'allozyme' in filetype:
			al.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')
		elif 'fresh' in filetype:
			fr.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')
		elif 'FFS' in filetype:
			ffs.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')
		elif 'rad' in filetype:
			rad.append(map)
			fh_out.write(s+'\t'+str(round(map,2))+'\n')

ave_fr = str(round(numpy.average(fr),2))
ave_al = str(round(numpy.average(al),2))
ave_ffs = str(round(numpy.average(ffs),2))
ave_rad = str(round(numpy.average(rad),2))

sd_fr = str(round(numpy.std(fr),2))
sd_al = str(round(numpy.std(al),2))
sd_ffs = str(round(numpy.std(ffs),2))
sd_rad = str(round(numpy.std(rad),2))

fh_out.write('\n')
fh_out.write('ave_mapping_fr'+'\t'+ave_fr+'\t'+ sd_fr +'\n')
fh_out.write('ave_mapping_al'+'\t'+ave_al+'\t'+ sd_al +'\n')
fh_out.write('ave_mapping_ffs'+'\t'+ave_ffs+'\t'+ sd_ffs +'\n')
fh_out.write('ave_mapping_rad'+'\t'+ave_rad+'\t'+ sd_rad +'\n')

