#!/usr/bin/env python -i

## Jazz18Heap.py version AA 01262020
## This program is a fast metasearch that will look for instances of sequence reads matching
## a reference in at least one k-mer.
## Jazz18Heap is intended to look for evidence of matches to a relatively short reference sequence
##    e.g., less than 100KB, but probably workable up to several MB in a large number of high throughput sequencing experiments
## Jazz18Heap is intended for finding relatively rare sequence (not common ones)
## Jazz18Heap doesn't substitute for many tools to align and track coverage.  It's main goal is rapid identification of potentially homologous sequences (with homology defined as a perfect match to a long sequence [e.g. 32bp for a 32b k-mer)
## Inputs are as follows (command line, Key=Value syntax)
##     ReferenceFile = <FastA file with list of sequences to match k-mers from>
##     ExcludeFile = <FastA file with sequences that will be excluded from matches>
##     DataFiles = <List of .fastq files, .fasta files, or NCBA-SRA accessions [e.g, SRR####]
##        Lists are comma delimited with no spaces in list
##           .fasta or .fastq files can be gzip compressed, although this somewhat slows the program down
##           .fasta and .fasta.gz files must have exactly one sequence per line (no multiline sequences)
##           For .fasta files to be downloaded from NCBI-SRA archive this entails the command line parameter --fasta 0
##        * Wildcards are allowed here as well, or list of files in a file with the extension .files
##        Providing a directory here will search all files in this directory or subdirectory for fasta and fastq data files
## Optional Parameters (willdefault to reasonable values if not set)
##     OutFileBase = <Character String to Label Output Files with>
##     ReportGranularity = <How many reads to process before reporting hit numbers (default is 1Million)>
##     SearchGranularity = <How much distance between k-mers to be examined in each data read
##         setting SearchGranularity=1 makes Jass18Heap look at every k-mer in every read
##         setting SearchGranularity=8 makes Jass18Heap look at every 8th k-mer in every read
##         higher numbers may miss a few hits but can greatly improve speed.
##         setting to a large number (999999) ensures only one k-mer will be looked up per read
##     SearchOffset = <Where to start jumping through each read for potential k-mers (zero means start at first base)
##     fastqdump = <Where to look for the fasterq-dump binary [program will look for this but if it is not found, the full path to fasterq-dump binary will need to be specified on the command line-- this is only important if you are downloading files from NCBI SRA as part of executing Jazz18Heap )
##     klen = <How long are the k-mers used>.  Default klen = 32 
##     snpAllow = <Set to True to allow a single mismatch in each k-mer (default is False)
##     Circular = <Set to True to force every Reference sequence to be treated as a circle>
##        Default : Uses the FastA name line-- if this line contains "Circular", the sequence is treated as a circle
## Multithreading: Jazz18Heap has the very primative multitasking ability to spawn a number
##     of derivative processes for a large number of data files to be scanned.  To use 16 Threads
##     Set MultiThreading=16 in the command line.
## Output:
##     Output consists of a log file with information on the run and a FastA file.
##     FastA files have information about each read in the ID line
##     ID Structure:
##          DataFileName;LineNumberInDataFile;ReferenceFileName;ReferenceSequenceName;NumberOfMatchedKmers
## Running the program:
##   Jazz18Heap only runs at full speed with a variant of Python (PyPY) that included a just-in-time compiler
##   Syntax Jazz18Heap RefFile=<MyRefFile> DataFiles=MyFastA1.fasta,MyFastA2.fasta.gz,MyFastQ*.fastq <Other_Parameters>

import gzip, os
from sys import argv, version
from time import time, strftime, localtime, asctime
from itertools import chain
from glob import glob
import subprocess
from random import choice


t0 = time()
now1 = strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())
Fn0 = []  ## List of reference files to search for sequences from
Gn0 = []  ## List of files containing sequences to exclude from k-mer match 
DataFile0 = []
OutFileBase0 = ''
klen1 = 0
ReportGranularity1 = 1000000  ## How often to report granularity, zero to provide less info in reports
SearchGranularity1 = 8 ## Granularity of search-- setting this to 1 searches ever k-mer in data files, setting to 10 will search every 10th k-mer
SearchOffset1 = 0 ## Granularity of search-- setting this to 1 searches ever k-mer in data files, setting to 10 will search every 10th k-mer
OneSnp1 = False ## Allow a single snp in each k-mer?
Circular1 = True ## Assume every sequence is circular?

MyCode1 = argv[0]
ai1 = 1

FastQDumpProgram1 = 'fasterq-dump' ## Location of a program that will download files from SRA if needed, can be reset by setting FastQDump=<path to fasterq-dump executable>.

MultiTask1 = False
BlockMultiCall1 = False
MultiCall1 = False
Multi1 = 0
Multi2 = 1

while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('ref') or a11.startswith('include'):
        Fn0.extend(a22.split(','))
    if a11.startswith('exclude'):
        Gn0.extend(a22.split(','))
    elif a11.startswith('data'):
        for a222 in a22.split(','):
            if os.path.isdir(a222):
                for D1,S1,FList1 in os.walk(a222):
                    for fc1 in FList1:
                        if fc1.endswith('.fastq') or fc1.endswith('.fastq.gz') or fc1.endswith('.fastq') or fc1.endswith('.fasta.gz'):
                            DataFile0.append(os.path.join(D1,fc1))
            elif '*' in a222:
                DataFile0.extend(glob(a222))
            elif a222.endswith('.files') or a222.endswith('.files.txt'):
                a2222 = open(a222,mode='r').read().split()
                for fc1 in a2222:
                    fc1 = fc1.strip()
                    if not(fc1):
                        continue
                    elif '*' in fc1:
                        DataFile0.extend(glob(fc1))
                    elif os.path.isdir(a222):
                        for D1,S1,FList1 in os.walk(a222):
                            for fc11 in FList1:
                                if fc11.endswith('.fastq') or fc11.endswith('.fastq.gz') or fc11.endswith('.fastq') or fc11.endswith('.fasta.gz'):
                                    DataFile0.append(os.path.join(D1,fc11))                    
                    else:
                        DataFile0.append(fc1)
            else:
                DataFile0.append(a222)
    elif a11.startswith('out'):
        OutFileBase0 = a22
    elif a11.startswith('k'):
        klen1 = int(a22)
    elif a11.startswith('reportgran'):
        ReportGranularity1 = int(a22)
    elif a11.startswith('searchgran'):
        SearchGranularity1 = int(a22)
    elif a11.startswith('searchoff'):
        SearchOffset1 = int(a22)
    elif a11.startswith('cir') :
        if a22.lower().startswith('f'):
            Circular1 = False
        else:
            Circular1 = True
    elif a11.startswith('snp') :
        if a22.lower().startswith('f'):
            OneSnp1 = False
        else:
            OneSnp1 = True
    elif a11.startswith('multi'):
        Multi1 = 0
        Multi2 = int(a22)
        MultiCall1 = True
    elif a11.startswith('instance'): 
        if '/' in a22:
            Multi1 = int(a22.split('/')[0])-1
            Multi2 = int(a22.split('/')[1])
            MultiTask1 = True
            MultiCall1 = False
            BlockMultiCall1 = True
    elif a11.startswith('fastqdump') or a11.startswith('fastq-dump') or a11.startswith('fasterqdump') or a11.startswith('fasterq-dump'):
        if os.path.isfile(a22) and not(os.path.isdir(a22)):
            FastQDumpProgram1 = a22
        elif os.path.isdir(a22):
            FastQDumpProgram1 = os.path.join(a22,'fasterq-dump')
            if not(os.path.isfile(FastQDumpProgram1)):
                FastQDumpProgram1 = os.path.join(a22,'bin','fasterq-dump')
                if not(os.path.isfile(FastQDumpProgram1)):
                    os.environ["PATH"]=os.getenv("PATH")+':'+a22
MultiMnemonic1 = 'J18H '
if MultiTask1 or MultiCall1:
    MultiMnemonic1 = 'J18H '+str(Multi1+1)+'/'+str(Multi2)+' '
if MultiCall1 and not(BlockMultiCall1):
    for j in range(1,Multi2):
        subprocess.Popen(['pypy',]+argv+['instance='+str(j+1)+'/'+str(Multi2)])
        
                                 
##Default Values for Debugging
if not Fn0:
    Fn0 = ['WuhanCorona0120PlusPolinton.fa',]
if not(DataFile0):
    DataFile0 = ['/Users/firelab08/Desktop/SRR3083808.fasta.gz',]
if not(OutFileBase0):
    OutFileBase0 = 'Jazz18Heap_I'+'{:02d}'.format(Multi1)+'_'+now1
if not(klen1):
    klen1 = 32


LogFile1="LogFile_"+OutFileBase0+'.tdv'
def LogNote1(note):
    LogFile=open(LogFile1,mode='a')
    LogFile.write(note+'\t'+'; Time='+"{0:.2f}".format(time()-t0)+' sec'+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' \r')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').strip(',') + '; Time='+"{0:.2f}".format(time()-t0)+' sec')
##Default Values for Debugging
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
LogNote1("Running Jazz18Heap with parameters:"+' '.join(argv)+' #Python Version'+version)
LogNote1("Python Flavor/Version: "+version)

if not('pypy' in version.lower()):
    LogNote1("Note that you should run Jazz18Heap with PyPy ('www.pypy.org') for full speed (otherwise speeds may be 10x sloer)") 
def LogOpeningFile1(FileID):
    LogNote1("Opening "+FileInfo1(FileID))
def LogClosingFile1(FileID):
    LogNote1("Closed "+FileInfo1(FileID))
def LogRunningFile1(FileID):
    LogNote1("Running "+FileInfo1(FileID))
LogRunningFile1(MyCode1)
if not Fn0:
    Fn0 = ['WuhanCorona0120.fa',]
    LogNote1('Reference File List set as debugging default to: '+str(Fn0)+ '. <RefFile= option on command line>')
if not(DataFile0):
    DataFile0 = ['/Users/firelab08/Desktop/SRR3083808.fasta.gz',]
    LogNote1('Data File List set as debugging default to: '+str(DataFile0)+ '. <DataFile= option on command line>')
if not(OutFileBase0):
    OutFileBase0 = 'CaughtReads_'+now1
    LogNote1('Output File Name Base set as debugging default to: '+str(OutFileBase0)+ '. <OutFile= option on command line>')
if not(klen1):
    klen1 = 32
    LogNote1('k-mer length set as debugging default to: '+str(klen1)+ '. <Klen= option on command line>')


D0 = {}  ## keys are k-mer sequences, values are tuples
         ## 0: Ordinal number of sequence that is hit
         ## 1: Ordinal position in sequence (one based +1 for sense, - for antisense)
         ## 2: Number of mismatches (presently 0 or 1)
NameList1 = []
FileList1 = []
SeqList1 = []
LenList1 = []
AntiList1 = []
snpD1 = {'G':'ATC','A':'TCG','T':'CGA','C':'GAT'}
def antisense1(seq):
    return seq.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
for nF0 in Fn0:
    Fx0 = open(nF0,mode='rU')
    LogOpeningFile1(Fx0)
    F0 = chain(Fx0,['>',])
    s1 = []
    n1 = os.path.basename(Fx0.name).split('.')[0]
    iN0 = 0
    for L0 in F0:
        if L0[0]=='>':           
            S1 = ''.join(s1).upper().replace('N','')  ##for now delete n's.  this willalso change some base numbers
            if S1:
                NameList1.append(n1)
                FileList1.append(nF0)
            if len(S1)>klen1:
                if ("circular" in n1.lower()) or Circular1:
                    S1 = S1+S1[:klen1-1]  #assume everything is a circle
                lD1 = len(S1)
                A1 = antisense1(S1)
                N1 = len(NameList1)-1
                LenList1.append(lD1)
                SeqList1.append(S1)
                AntiList1.append(A1)                
                for i in xrange(len(S1)-klen1+1):
                    K1 = S1[i:i+klen1]
                    if (K1 in D0) and D0[K1][2]==0:
                        continue
                    D0[K1] = (N1,i+1,0)
                    D0[antisense1(K1)] = (N1,-i-1,0) ## antisense of a k-mer at position p will have the position -p.  Note all positions are one-based 
                    if OneSnp1:
                        for j in xrange(klen1):
                            for k in xrange(3):
                                mutK1 = K1[:j]+snpD1[K1[j]][k]+K1[j+1:]
                                if not(mutK1 in D0):
                                    D0[mutK1] = (N1,i+1,1)
                                    D0[antisense1(mutK1)] = (N1,-i-1,1)
            n1a = L0[1:].strip().split()
            n2a = []
            for n0a in n1a:
                if not("=" in n0a) or ("range=" in n0a):
                    n2a.append(n0a.replace("range=",""))
            n1 = '_'.join(n2a)
            s1 = []
        else:
            s1.append(L0.strip())                         
    Fx0.close()
    LogNote1('Added k-mer dictionary from '+nF0)
LogNote1('Finished making initial index')

for nF0 in Gn0:  ## currently set up to only remove perfect matches
    Fx0 = open(nF0,mode='rU')
    LogOpeningFile1(Fx0)
    F0 = chain(Fx0,['>',])
    s1 = []
    n1 = 'Default_'+nF0
    iN0 = 0
    for L0 in F0:
        if L0[0]=='>':           
            S1 = ''.join(s1)
            if len(S1)>klen1:
                if ("circular" in n1.lower()) or Circular1:
                    S1 = S1+S1[:klen1-1]  #assume everything is a circle
                A1 = S1.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()[::-1]
                for i in xrange(len(S1)-klen1+1):
                    if S1[i:i+klen1] in D0:
                        del(D0[S1[i:i+klen1]])
                    if A1[i:i+klen1] in D0:
                        del(D0[A1[i:i+klen1]])
        else:
            s1.append(L0.strip())                         
    Fx0.close()
    LogNote1('Masked from k-mers from '+nF0)
LogNote1('Finished filtering index')

TempFileUID1 = ''
for i in range(12):
    TempFileUID1 += choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')
ScratchPath1 = os.path.join(os.getcwd(), 'Jazz18HeapTempFiles'+TempFileUID1)

DataFile1 = []
for df0 in DataFile0[Multi1::Multi2]:
    if os.path.isfile(df0) and not(df0.endswith('.sra')):
        DataFile1.append(df0)
    else:
        LogNote1(df0+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
        LogNote1("Preparing to download sequence read set "+df0+" from NCBI")
        try:
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'version'])
        except:
            LogNote1("Searching for a version of fasterq-dump that will run; if this fails, you may need to redownload the program and unzip the archive, also add FastQDumpProgram=<path to program> to command line")
            os.environ["PATH"]=os.getenv("PATH")+':./:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/Applications:~/Downloads'
            for fqd1 in os.getenv("PATH").split(':'):
                if os.path.isdir(fqd1):
                    for fqd2 in os.listdir(fqd1):
                        if fqd2.startswith('sratoolkit'):
                            fqd3 = os.path.join(fqd1,fqd2)
                            if os.path.isdir(fqd3):
                                fqd4 = os.path.join(fqd3,'bin','fasterq-dump')
                                if os.path.isfile(fqd4):
                                    if versioner1(fqd2) > versioner1(fqdVersion1):
                                        fqdVersion1 = fqd2
                                        FastQDumpProgram1 = fqd4
                                        TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'version'])
        LogNote1("Trying presumed fasterq-dump program file located at "+FastQDumpProgram1)
        if not(os.path.isdir(ScratchPath1)):
            os.mkdir(ScratchPath1)
        TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                 '--fasta',
                                                 '0',
                                                 '--origfmt',
                                                 '--outdir',
                                                 ScratchPath1,
                                                 df0])
        LogNote1("Result of "+df0+" NCBI Download " + TryFastQDump1)
        PresumptiveFilePath1 = os.path.join(ScratchPath1,df0+'.fasta')
        if os.path.isfile(PresumptiveFilePath1):
            DataFile1.append(PresumptiveFilePath1)
        else:
            LogNote1('Looks Like fasterq-dump failed for '+df0)

OutFile1 = open('CaughtReads_'+OutFileBase0+'.fa',mode='w') ## FastA File
F2 = open('FileByFile_'+OutFileBase0+'.tdv',mode='w') ## FastA File
F2.write('File\tTotalReads\tTotalBases\tHitReads\tHitKmers\tMatchList\n')
OrientA1 = ['b','s','a']
LogOpeningFile1(OutFile1)
TotalReads1 = 0
TotalHitReads1 = 0
TotalBases1 = 0
TotalHitKmers1 = 0
for iD0,nD0 in enumerate(DataFile1):
    Roton1 = 1
    HotLine1 = 2
    if nD0.lower().endswith('fastq') or nD0.lower().endswith('fastq.gz'):
        Roton1 = 4
        HotLine1 = 2
    if nD0.lower().endswith('fasta') or nD0.lower().endswith('fasta.gz'):
        Roton1 = 2
        HotLine1 = 2
    if nD0.endswith('gz'):
        if version.startswith('2.'):
            F1 = gzip.open(nD0,mode='r')
            
        else:
            F1 = gzip.open(nD0,mode='rt')
    else:
        if version.startswith('2.'):
            F1 = open(nD0,mode='rU')
        else:
            F1 = open(nD0,mode='r')
    if ReportGranularity1:
        LogOpeningFile1(nD0)
    F1r = F1 ## may be faster with .read().splitlines()
    ReportInterval1 = ReportGranularity1 * Roton1
    FileReads1 = 0
    FileHitReads1 = 0
    FileBases1 = 0
    FileHitKmers1 = 0
    RotonCounter1 = 0
    LineCounter1 = 0
    HitD1 = {}
    ReadFileMnemonic1 = os.path.basename(F1.name).split('.')[0]
    nL1=-1
    for L1 in F1r:
        nL1+=1
        RotonCounter1 += 1
        if RotonCounter1 == HotLine1:
            RotonCounter1 = 0
            R1 = L1.strip()
            len1 = len(R1)
            m1 = 0
            for i in xrange(SearchOffset1,len1-klen1+1,SearchGranularity1):
                q1 = R1[i:i+klen1]
                if q1 in D0:
                    if m1 == 0:
                        h1 = D0[q1]
                        h = h1[0]
                        m1 = 1
                        p1 = h1[1]
                        if p1>0:
                            p1 -= i
                        else:
                            p1 = -(-p1+klen1+i)
                    elif D0[q1][0]==h:
                        m1 += 1
            if m1>0:
                FileHitReads1 += 1
                FileHitKmers1 += m1
                if p1>0:
                    OutFile1.write('>'+
                               ReadFileMnemonic1+
                               '_Read='+str(nL1//Roton1)+
                               '_Source='+FileList1[h]+
                               '_Seq='+NameList1[h]+
                               '_Pos='+str(p1)+'s'+
                               '_KMatches='+str(m1)+
                               '\r')
                else:
                    OutFile1.write('>'+
                               ReadFileMnemonic1+
                               '_Read='+str(nL1//Roton1)+
                               '_Source='+FileList1[h]+
                               '_Seq='+NameList1[h]+
                               '_Pos='+str(p1)+'a'+
                               '_KMatches='+str(m1)+
                               '\r')
                OutFile1.write(R1+'\r')
                if not h in HitD1:
                    HitD1[h]={}
                if p1>0:
                    if not((p1,m1,1) in HitD1[h]):
                        HitD1[h][(p1,m1,1)] = 0
                    HitD1[h][(p1,m1,1)] += 1
                else:
                    if not((-p1,m1,-1) in HitD1[h1[0]]):
                        HitD1[h][(-p1,m1,-1)] = 0
                    HitD1[h][(-p1,m1,-1)] += 1
            FileReads1 += 1
            FileBases1 += len1
        LineCounter1 += 1
        if ReportGranularity1 and LineCounter1==ReportInterval1:
            LineCounter1 = 0
            LogNote1(MultiMnemonic1+'Completed Read '+
                     str(1+nL1//Roton1)+
                     ' from '+
                     ReadFileMnemonic1 +
                     '  Bases='+
                     str(FileBases1)+
                     '  HitReads='+
                     str(FileHitReads1)+
                     '  HitKmers='+
                     str(FileHitKmers1))                     
    if ReportGranularity1:
        LogClosingFile1(F1)
    F1.close()
    rdT1 = ''
    for i00 in sorted(list(HitD1.keys())):
        rdT1+=NameList1[i00]+':'
        for (p1,m1,o1) in sorted(list(HitD1[i00].keys())):
            rdT1 += 'p'+str(p1)+OrientA1[o1]+'_m'+str(m1)
            if HitD1[i00][(p1,m1,o1)]>1:
                rdT1 += '('+str(HitD1[i00][(p1,m1,o1)])+')'
            rdT1 += ' '
    rdT1 = rdT1.strip()
    rdT2 =''
    if rdT1:
        rdT2 = ' Matches='+rdT1
    LogNote1(MultiMnemonic1+'Finishing File '+str(iD0)+': '
             +ReadFileMnemonic1+
             ' Reads='+
             str(FileReads1)+
             ' Bases='+
             str(FileBases1)+
             ' HitReads='+
             str(FileHitReads1)+
             ' HitKmers='+
             str(FileHitKmers1)+
             rdT2)
    if ReportGranularity1:
        LogNote1(MultiMnemonic1+'Project Progress: '
         +'AllFilesSoFar'+
         ' Reads='+
         str(TotalReads1)+
         ' Bases='+
         str(TotalBases1)+
         ' HitReads='+
         str(TotalHitReads1)+
         ' HitKmers='+
         str(TotalHitKmers1))
    F2.write(nD0+
             '\t'+
             str(FileReads1)+
             '\t'+
             str(FileBases1)+
             '\t'+
             str(FileHitReads1)+
             '\t'+
             str(FileHitKmers1)+
             '\t'+
             rdT1+'\r')
    TotalReads1 += FileReads1
    TotalBases1 += FileBases1
    TotalHitReads1 += FileHitReads1
    TotalHitKmers1 += FileHitKmers1
    if nD0.startswith(ScratchPath1) and nD0.endswith('.fasta'):
        os.remove(nD0)
    
LogNote1(MultiMnemonic1+'Finishing Run: '
         +'AllFiles'+
         ' Reads='+
         str(TotalReads1)+
         ' Bases='+
         str(TotalBases1)+
         ' HitReads='+
         str(TotalHitReads1)+
         ' HitKmers='+
         str(TotalHitKmers1))
F2.write('AllFiles'+
             '\t'+
             str(TotalReads1)+
             '\t'+
             str(TotalBases1)+
             '\t'+
             str(TotalHitReads1)+
             '\t'+
             str(TotalHitKmers1)+'\r')
LogClosingFile1(F2)
F2.close()
LogClosingFile1(OutFile1)
OutFile1.close()        
LogNote1('Finished running '+' '.join(argv))
try:
    LogNote1(open(MyCode1,mode='r').read())
except:
    pass

##  Jazz18Heap a fast Python Script for finding matches to a probe sequence in large numbers of
##  High throughput sequencing output files
##  Copywrite 2020, Andrew Fire and Stanford University
##  With thanks to Dae-Eun Jeong, Loren Hansen, Matt McCoy, Nimit Jain, Massa Shoura, Karen Artiles, Lamia Wahba, Nelson Hall
     
    
        


        
        
        
