import os
import shutil

'''
zip cleaned and unfiltered dirs

'''
backup = 'cleaned_unfiltered_backup'
if not os.path.exists(backup):
	os.mkdir(backup)
files = []
for dir in os.listdir('.'):
	if os.path.isdir(dir):
		if not dir.startswith('cleaned') and not dir.startswith('gyro'):
			print 'moving ', dir
			files.append(dir)
			#shutil.make_archive(dir, 'gztar', dir)
			shutil.move(dir,backup)

#print 'zipping ', backup
#shutil.make_archive(backup, 'gztar', backup)

for f in files:
	print 'making dir ', f
	if not os.path.exists(f):
		os.mkdir(f)

for f in files:
	print 'moving and renaming files ', f
	forward = f+'_clean-READ1.tagged_filter.fastq'
	forward_single = f+'_clean-READ1-single.tagged_filter.fastq'
	reverse = f+'_clean-READ2.tagged_filter.fastq'
	reverse_single = f+'_clean-READ2-single.tagged_filter.fastq'
	if os.path.exists(forward):
		shutil.move(f+'_clean-READ1.tagged_filter.fastq',f+'/'+f+'_clean-READ1.fastq')
	if os.path.exists(reverse):
		shutil.move(f+'_clean-READ2.tagged_filter.fastq',f+'/'+f+'_clean-READ2.fastq')
	if os.path.exists(forward_single):
		shutil.move(f+'_clean-READ1-single.tagged_filter.fastq',f+'/'+f+'_clean-READ1-single.fastq')
	if os.path.exists(reverse_single):
		shutil.move(f+'_clean-READ2-single.tagged_filter.fastq',f+'/'+f+'_clean-READ2-single.fastq')

for f in files:
	for i in os.listdir('.'):
		if i.startswith(f):
			if os.path.isdir(i):
				pass
			else:
				os.remove(i)
