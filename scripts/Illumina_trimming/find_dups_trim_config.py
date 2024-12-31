i7 = 'CCGAGCCCACGAGAC*ATCTCGTATGCCGTCTTCTGCTTG'
i5 = 'AATGATACGGCGACCACCGAGATCTACAC*TCGTCGGCAGCGTC'

config = 'trimmomatic_config.txt'
fh_conf = open(config,'r')
fh_out = open('testing.txt','a')
samples = []
bar7 = []
bar5 = []
input_list = []
extra_list = []
input_list_filt = []
dup_list = []

for line in fh_conf:
	line = line.strip()
	if line.startswith('['):
		pass
	elif ':_' in line:
		samples.append(line)
	elif line.startswith('i7'):
		bar7.append(line)
	elif line.startswith('i5'):
		bar5.append(line)
for a,b,c in zip(samples,bar7,bar5):
	input_list.append(a+' ' + b + ' ' + c)

###create way of grabbing the dup samples so I can come back to them
for sample in input_list:
	s_name = sample.split(' ')[0].split(':')[0]
	if s_name not in extra_list:
		extra_list.append(s_name)
		input_list_filt.append(sample)
	else:
		dup_list.append(sample)
print 'input_list',len(input_list)
print 'filt', len(input_list_filt)
print 'dup', len(dup_list)
print 'diff', len(input_list)-len(input_list_filt)
for item in dup_list:
	fh_out.write(item+'\n')