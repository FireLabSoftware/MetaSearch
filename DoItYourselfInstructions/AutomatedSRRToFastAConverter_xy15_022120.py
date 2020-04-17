#!/usr/bin/env python -i

##
## -> AutomatedSRRToFastAConverter_xy15_022120
## -> AutomatedSRRToFastAConverter takes a directory with NCBI-SRA files
##      and converts all .sra files into fasta format using fastq-dump.
##      Humble Note: You can also do this with
##      a single command line argument ./fastq-dump /mydir/*.sra
##      but that may lose some flexibility in the call
##
## -> Program operates with multple threads at once and operates both from 
##      beginning to end and from end to beginning of a list of files
##      so a single failure won't as easily stop the entire ship
##
## -> Assumes that SRA tools is installed on your system and that
##      fastq-dump is available (fasterq-dump could also be used if you can get
##      the fasterq-dump dependencies to work.
##
## -> Tested with Python 2.7 (will maybe run with Python 3 but not tested)
##
## -> Call Format
##       python AutomatedSRRToFastAConverter_ab_021620 <parameters/options>
##           Required Parameters:
##               SourceDirectory = <Input directory>  *You'll need to provide*
##               DestinationDirectory = <Destination_Directory> *Default is to use current directory*
##           Options:
##               FastQDumpProgram = <Path_to_fastq-dump> *Default is to look in current directory*
##                 or FasterQDumpProgram = <Path_to_fasterq-dump> *Default is to look in current directory*
##               Threads = <Number of Parallel threads to try> *Default is 16*
## -> Copyright Fire Lab, Stanford University 021620, No Guarantees of any sort
##
## -> Version AC 021720 (First version with command line interface)
##
from sys import argv
import os
import subprocess
from time import time,strftime,localtime
from glob import glob

Threads1 = 16
Source1 = './*.sra'
FastQDumpProgram1 = './fastq-dump' ##'./fasterq-dump' ## Can use either program if it works on your system if the program works (fasterq-dump may have additional dependencies)
Destination1 = './'
t0 = time()
ThisProgram1 = argv[0]
Meta1 = True
ai1=1
mycore1 = 0
myorder1 = 'F'
while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if (a1.startswith('help') or ('-help' in a1)) and Meta1:
        mycode1 = open(argv[0]).read().splitlines()
        for h1 in mycode1:
            if h1.startswith('##'):
                print(h1[2:])
            elif h1.startswith('import'):
                break
        exit()
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('thread'):
        Threads1 = int(a22)
    elif a11.startswith('fast'):
        FastQDumpProgram1 = a22
    elif a11.startswith('destination'):
        Destination1 = a22
    elif a11.startswith('source'):
        Source1 = a22 
    elif a11.startswith('mycore'):
        mycore1 = int(a22)
        Meta1 = False
    elif a11.startswith('myorder'):
        myorder1 = a22
        Meta1 = False

now1 = strftime("%m_%d_%y_%H_%M_%S",localtime()) ## String Representation of Current Time

if os.path.isdir(Source1):
    DataIterator1 = os.walk(Source1)
else:
    DataIterator1 = glob(Source1)

A = []
cores1 = Threads1/2

if Meta1:
    for i in range(cores1):
        A.append(subprocess.Popen(['python',
                                   ThisProgram1,
                                   "Source="+Source1,
                                   "Destination="+Destination1,
                                   "fastqdump="+FastQDumpProgram1,
                                   "mycore="+str(i),
                                   "myorder=F"]))
    for i in range(cores1):
        A.append(subprocess.Popen(['python',
                                   ThisProgram1,
                                   "Source="+Source1,
                                   "Destination="+Destination1,
                                   "fastqdump="+FastQDumpProgram1,
                                   "mycore="+str(i),
                                   "myorder=R"]))
    for a in A:
        a.wait()
    exit()

else:
    MyID1 = str(mycore1)+myorder1
    q=0
    G = open('SraToFastAOutputInfo'+MyID1+'.txt',mode='w')
    s1 = "Starting " +ThisProgram1 + " in core mode "+MyID1+". TimeInSeconds=%.2f"%(time()-t0)
    G.write(s1+'\r'); print(s1)
    DataIterator1 = list(DataIterator1)
    if myorder1 == 'R':
        DataIterator1 = DataIterator1[::-1]
    for dsf1 in DataIterator1:
        if len(dsf1)==3:
            (d,s,f)=dsf1
        else:
            d=''
            f = (dsf1,)
        if myorder1 == 'R':
            f = f[::-1]
        for fn1 in f:
            if not(fn1.endswith('.sra')):
                continue
            q+=1
            if q % cores1 == mycore1:
                S = os.path.join( d , fn1 )
                if os.path.isfile(os.path.join(Destination1,fn1.split('.')[0]+'_1.fasta.gz')):
                    continue
                s1=MyID1," File#=",q,".  Dataset=",fn1,". TimeInSeconds=",'%.2f'%(time()-t0)
                s3=''.join(map(str,s1))
                G.write(s3+'\r'); print(s3)
                s2=subprocess.check_output([FastQDumpProgram1,
                                            '--split-files',
                                            '--fasta',
                                            '0',
                                            '--origfmt',
                                            '--gzip',
                                            '--outdir',
                                            Destination1,
                                            S])
                s3 = s2.decode("utf-8")
                G.write(s3+'\r'); print(s3)
G.close()
