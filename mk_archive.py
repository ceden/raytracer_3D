
import os,glob,tarfile

def always_true(file):
   if file == 'bin': return False
   return True

archive_name = 'tmp.tar.gz'
print 'creating archive ',archive_name 
tar = tarfile.open(archive_name,'w:gz')
tar.add('bin',exclude=always_true)
for file in glob.glob('bin/*.ipynb'): tar.add(file)
#tar.add('doc/doc1.pdf')
tar.add( 'for_config/Makefile' )
for file in glob.glob('for_config/*.f90'): tar.add(file)
tar.add( 'for_src/Makefile' )
for file in glob.glob('for_src/*.f90'): tar.add(file)
tar.add( 'make.inc' )

tar.close()

if 1:
  tar = tarfile.open(archive_name,'r:gz')
  tar.list()
  tar.close()

