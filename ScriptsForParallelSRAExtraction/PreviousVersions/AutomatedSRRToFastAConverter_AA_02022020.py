import sys, os, subprocess, time
cores1 = 8
SRADir1 = '/media/fl19/HDD8TB01/AMPData/sra/'
FastQDumpProgram1 = './fastq-dump' ##'./fasterq-dump'
DataIterator1 = os.walk(SRADir1)
FastADir1 = '/media/fl19/HDD8TB01/AMPData/fasta/'
t0 = time.time()
ThisProgram1 = sys.argv[0]
if len(sys.argv)==1:
    r1='Meta'
    o1='Bidirectional'
else:
    r1=int(sys.argv[1])
    o1=sys.argv[2]
if r1=='Meta':
    for i in range(cores1):
        subprocess.Popen(['python',ThisProgram1,str(i),'F'])
    for i in range(cores1):
        subprocess.Popen(['python',ThisProgram1,str(i),'R'])
else:
    MyID1 = str(r1)+o1
    q=0
    G = open('SraToFastAOutputInfo'+str(r1)+o1+'.txt',mode='w')
    s1 = "Starting " +ThisProgram1 + " in core mode "+str(r1)+", orientation "+o1+". TimeInSeconds=%.2f"%(time.time()-t0)
    G.write(s1+'\r'); print(s1)
    for d,s,f in DataIterator1:
        if o1 == 'R':
            f = f[::-1]
        for fn1 in f:
            q+=1
            if q % cores1 == r1:
                S = os.path.join( d , fn1 )
                if os.path.isfile(os.path.join(FastADir1,fn1.split('.')[0]+'_1.fasta.gz')):
                    continue
                if os.path.isfile(os.path.join(FastADir1,fn1.split('.')[0]+'_2.fasta.gz')):
                    continue
                s1=MyID1," File#=",q,".  Dataset=",fn1,". TimeInSeconds=",'%.2f'%(time.time()-t0)
                s3=''.join(map(str,s1))
                G.write(s3+'\r'); print(s3)
                s2=subprocess.check_output([FastQDumpProgram1,
                                            '--split-files',
                                            '--fasta',
                                            '0',
                                            '--origfmt',
                                            '--gzip',
                                            '--outdir',
                                            FastADir1,
                                            S])
                s3 = s2.decode("utf-8")
                G.write(s3+'\r'); print(s3)
    s4 = "#####Finishing " +ThisProgram1 + " in core mode "+str(r1)+", orientation "+o1+". TimeInSeconds=%.2f"%(time.time()-t0)
    G.write(s4+'\r'); print(s4)
    G.close()
