#!/usr/bin/env python -i

##
## -> AllInOnePrefetcher_xy15_022120
## -> AllInOnePrefetcher takes a list of NCBI-SRA accession numbers
##      and downloads these to a specific folder
## -> Program operates with multple threads at once and tries continuously
##      to download each requested file to avoid download failures
## -> Assumes that SRA tools and Aspera Connect are installed on your system
## -> Some files are protected by IRB and the program will generally skip these.
##      To try to download, set irb=True (not tested)
## -> Very large files may be unwanted.  Current Maximum download is set to 80G
##      (larger files will be skipped),  To change this, set MaxSize=##G
## -> Tested with Python 2.7 (will maybe run with Python 3 but not tested)
## -> Call Format
##       python AllInOnePrefetcher <parameters/options>
##         Required Parameters:
##           ListFile=<SRARunTable_file>   *You'll need to provide a list*
##           ListFile can be in any format as long as the first item is an SRA accession##
##         Options:
##           PrefetchProgram = <Path_to_Prefetch> *Default is to look in current directory*
##           DestinationDirectory = <Destination_Directory> *Default is to use sra-tools setting*
##           Threads = <Number of Parallel threads to try> *Default is 16*
##           MaxSize = <Maximum download size in Gigabytes> *Default is 80*
##           IRB = <> *Attempt to download IRB-protected files? [default is False]*
##           Delimiter = <> Column delimiter used in ListFile (generally tab or comma)
##               Should autodetect for recent NCBI SRARunTable formats
##               Otherwise can explicitly set with Delimiter=tab or Delimiter=comma
##           RunColumn = <> Column # for RunID in ListFile.
##               Should autodetect if you have retained the header row in the SRATable
##               Defaults 0 comma-delimited (new RunTable fmt), 32 tab-delimited (old fmt)
##           ExclusionFile = <> File with a list of archive names to exclude (e.g. already analyzed such as a file-by-file list)
##           ExclusionDir =  <> A directory of SRA archive files to skip 
## -> Copyright Fire Lab, Stanford University 02-2020, No Guarantees of any sort
##
from sys import argv
import os
import subprocess
from time import time,strftime,localtime

Meta1 = True ## As default, running the program spawns child process (for a total of Threads1 processes).  This variable is true for the main running of the program
Threads1 = 16  ## Number of threads to use in the system (up to twice number of cores may be useful in some circumstances
PrefetchProgram1 = './prefetch' ## Default location for prefetch program-- can be set in command line to any location
TryToDownloadIRBFiles1 = False ## Set to False unless you have authorization to download datasets with '-IRB-' in metadata
MaxSize1 = '80G' ## Maximum Size of file to attempt download  
Destination1 = './' ## Destination for .sra files  
ExclusionD1 = {}

def FullSplit1(s0):
	s1 = s0.replace('.',' ').replace(',',' ').replace('/',' ').replace(';',' ').replace(':',' ').replace('_',' ').split()
	s2 = [w2 for w2 in s1 if w2]
	return s2
	
## Get input from users as to what they want to obtain and where they want to put it (and how many processor threads to use)
ai1=1
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
    elif a11.startswith('prefetch'):
        PrefetchProgram1 = a22
    elif a11.startswith('delimiter'):
        if a22.lower().startswith('t') or a22.lower().startswith('\t') or a22.lower().startswith('\\t'):
            Delimiter1 = '\t'
        elif a22.lower().startswith('c') or a22.lower().startswith(',') or a22.lower().startswith('\\c'):
            Delimiter1 = '\c'
        else:
            Delimiter1 = a22
    elif a11.startswith('runcolumn'):
        RunColumn1 = int(a22)        
    elif a11.startswith('destination'):
        Destination1 = a22
    elif a11.startswith('list') or a11.startswith('source') or a11.startswith('data'):
        ListFile1 = a22
    elif a11.startswith('mycore'):
        mycore1 = int(a22)
        Meta1 = False
    elif a11.startswith('exclu'):
		for a222 in a22.split(','):
			if os.path.isdir(a222):
				for (d,s,f) in os.walk(a222):
					for fn1 in f:
						if fn1.endswith('.sra'):
							ExclusionD1[fn1.split('.sra')[0]] = 0
			elif os.path.isfile(a222):
				if a222.endswith('.sra'):
					ExclusionD1[a222.split('.sra')[0]] = 0
				else:
					for L1 in open(a222,mode='rU'):
						for w1 in FullSplit1(L1):
							ExclusionD1[w1] = 0
			else:
				for fx1 in glob(a222):
					for fx2 in FullSplit1(fx1):
						ExclusionD1[fx2] = 0
    elif a11.startswith('maxsize'):
        try:    
            MaxSize1 = str(int(a22))+'G'
        except:
            MaxSize1 = a22
    elif a11.startswith('irb') :
        if a22.lower().startswith('f'):
            TryToDownloadIRBFiles1 = False
        else:
            TryToDownloadIRBFiles1 = True


now1 = strftime("%m_%d_%y_%H_%M_%S",localtime()) ## String Representation of Current Time
A = []  ## Keep the subprocess objects in an array for potential debugging

if Meta1:
    for i in range(Threads1):
        print("Starting Thread: "+str(i))
        A.append(subprocess.Popen(('python',)+tuple(argv)+('mycore='+str(i),)))
    for a in A:
        a.wait()
    exit()
r = int(mycore1)
t0 = time()
r1 = str(r)
if len(r1)==1:
    r1 = '0'+r1
F = open(ListFile1,mode='rU').read().splitlines()
G = open('SraOutputInfo_'+r1+'_'+now1+'.txt',mode='w')
q=0
Caught1 = 1
Iteration1 = 1
Delimiter1 = 'x'
while Caught1>0:
    Caught1 = 0
    print len(F)
    for L in F:
        if Delimiter1=='x':
            if L.lower().startswith('run\t'):
                RunColumn1 = 0
                Delimiter1 = '\t'
                continue
            elif L.lower().startswith('run,'):
                RunColumn1 = 0
                Delimiter1 = ','
                continue
            elif ('\t' in L):
                Delimiter1 = '\t'
                V = L.lower().split('\t')
                if 'run' in V:
                    RunColumn1 = V.index('run')
                    continue
                else:
                    RunColumn1 = 32
            else:
                Delimiter1 = ','
                V = L.lower().split(',')
                if 'run' in V:
                    RunColumn1 = V.index('run')
                    continue
                else:
                    RunColumn1 = 0
        if '-IRB-' in L and not(TryToDownloadIRBFiles1):
            continue
        q+=1
        if q%Threads1==r:
            S = L.split(Delimiter1)[RunColumn1].strip()  ## Assume the SRA accession is the first field in each line
            if os.path.isfile(os.path.join(Destination1,S+'.sra')) or (S in ExclusionD1) or (S+'.sra' in ExclusionD1):
                continue
            Caught1+=1
            s1="Start File Number=",q,".  Dataset=",S,". TimeInSeconds=",'%.2f'%(time()-t0)
            s3=''.join(map(str,s1))
            G.write(s3)
            print(s3)
            s2 = ''
            try:
                s2=subprocess.check_output([PrefetchProgram1,'--output-directory',Destination1,'--max-size',MaxSize1,S]).decode("utf-8")
            except:
                pass
            G.write(s2+'\r')
            print(s2 +'\r')
    print("Finished Iteration *"+str(Iteration1)+"* This run caught "+str(Caught1)+" files from this set that had previously not been downloaded.")
    Iteration1 += 1
G.close()
