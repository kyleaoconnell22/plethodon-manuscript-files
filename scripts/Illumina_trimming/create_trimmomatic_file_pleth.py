'''
these three were not found
('Sample', 'cylindraceus-IGG-Fresh-RCB103', 'not found in control file. Skipped.')
('Sample', 'yonahlossee-SG-Fresh-', 'not found in control file. Skipped.')
('Sample', 'cylindraceus-IGG-Fresh-RCB102', 'not found in control file. Skipped.')

#Also, it can't handle the - so that needs to be removed
'''
fh1 = open('sample_list_all.txt','r')
fh3 = open('raw_samples_pleth.txt','r')
fh_out = open('trimmomatic_config.txt','a')
barcodes = []
raw_samples = []
test = []
for line in fh1: 
	line = line.strip()
	barcodes.append(line)
for line in fh3:
	line = line.strip()
	raw_samples.append(line)

fh_out.write('[adapters]'+'\n'+'i7:CCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG' + \
	'\n'+'i5:AATGATACGGCGACCACCGAGATCTACAC*TCGTCGGCAGCGTC'+'\n'+'\n')

fh_out.write('[names]'+'\n')

for item in raw_samples:
	for b in barcodes:
		#print item.split('-')[3], b.split('\t')[0]
		if item.split('_')[3] == b.split('\t')[0]:
			fh_out.write('{0}:_\n'.format(item))

fh_out.write('\n'+'[barcodes]'+'\n')


for item in raw_samples:
	for b in barcodes:
		#print item.split('-')[3], b.split('\t')[0]
		if item.split('_')[3] == b.split('\t')[0]:
			barcode = b.split('\t')[1]
			fh_out.write('i7-{0}:{1}\n'.format(item,barcode.split('-')[0]))
			fh_out.write('i5-{0}:{1}\n'.format(item,barcode.split('-')[1]))