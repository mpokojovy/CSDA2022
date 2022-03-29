import subprocess

f = open('jobsubmit','r')
js = f.read()
f.close()

for x in range(1,2):
   fo = open('jobsubmit_chunk_%s'%str(x),'w')
   fo.write(js)
   fo.write('\nR CMD BATCH --quiet --vanilla --no-restore chunk.%s.R > chunk.%s.out\n'%(str(x),str(x)))
   fo.close()

   process = subprocess.Popen('sbatch jobsubmit_chunk_%s'%str(x), stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
   process.wait()

   print(process.stdout.readlines())
   print(process.stderr.readlines())


