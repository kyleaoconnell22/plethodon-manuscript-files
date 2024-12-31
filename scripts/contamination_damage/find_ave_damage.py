import numpy

fh5=open('5pCtoT_freq.txt','r')
fh3=open('3pGtoA_freq.txt','r')

a = []
b = []

for line in fh5:
	line=line.strip()
	if not line.startswith('pos'):
		a.append(float(line.split('\t')[1]))

for line in fh3:
	line=line.strip()
	if not line.startswith('pos'):
		b.append(float(line.split('\t')[1]))

print 'ave 3p damage = ', numpy.mean(b)
print 'sd 3p damage = ', numpy.std(b)
print 'ave 5p damage = ', numpy.mean(a)
print 'sd 5p damage = ', numpy.std(a)

