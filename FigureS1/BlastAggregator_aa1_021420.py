#!/usr/bin/env python -i

## Version BlastAggregator_aa1_021420

## Blast Aggregator: Input is a blast search result from a short read archive probed with a reference sequence, along with reference sequence
## Output are FastA files with an approximate version of the input file with all regions that are covered in the blast search set to the SRA read sequences
## 'Democratic' resolution of conflicts where two reads disagree (with the winner in lower case if there is any disagreement
## Takes a standard BlastN alignment output from the web as of 02/20 as input, along with a FastA file of the reference used to generate the Blast Alignments
## Democratic resolutions use several factors below for tiebreaking
## Two output files, one containing insertions from the Blast alignments (so numbering may be different) and one without insertions
## Copyright Fire Lab, Stanford University 021020, No Guarantees of any sort

## Call Format python BlastAggregator##.py RefFile AlignFile

from sys import argv

## Input Reference File
RefFile = argv[1]  ## Original Reference File for the search in FastA format (one sequence only please)
AlignFile = argv[2] ## File of alignments from NCBI blast output on web

## Input Options
SwapQuerySubject1 = False ## True to swap the query and subject designations

## Parameters to resolve "tie votes" where two different sequences are present in equal numbers
ReferenceBias1 = 1.0/3  ## Upvote possibilities that match the reference
MultalignBias1 = -1.0/7 ## Downvote reads that map to more than one place
StartBias1 = [-1/11.0, -1/13.0, -1/17.0, -1/19.0, -1/23.0]  ##  slight downvote to first few in a read   
EndBias1 = [-1/43.0, -1/41.0, -1/37.0, -1/31.0, -1/29.0]  ##  slight downvote to last few in a read   
PosBiasDivisor1 = 104729.0  ## Very slight downvote to later bases in a read   

## Minor output formatting
LineLengthOut1 = 100
ReplaceDsWithDashes1 = True
LineEnd1 = '\r' ## Choose your line endings

def PosBias1(Qp):  
    return -Qp/PosBiasDivisor1

for x in argv[3:]:
    if x.lower().startswith('swap'):
        SwapQuerySubject1 = True

ASB11="AaCcN*NgGtT"
filterplusupper=''
for i in range(256):
    if chr(i) in ASB11:
        filterplusupper+=chr(i).upper()
    else:
        filterplusupper+=' '
filterplus=''
for i in range(256):
    if chr(i) in ASB11:
        filterplus+=chr(i)
    else:
        filterplus+=' '

def stripseq(s):
    '''return a filtered version of any sequence)'''
    return s.translate(filterplus).replace(' ','')
def stripsequpper(s):
    '''return a filtered version of any sequence)'''
    return s.translate(filterplusupper).replace(' ','')

BaseD1 = {'G':0, 'A':1, 'T':2, 'C':3, 'D':4, '-':4, 'g':0, 'a':1, 't':2, 'c':3, 'd':4}

ref = open(RefFile,mode='rU').read()
a = open(AlignFile,mode='rU').read()
rs1 = ''
rn1 = ''
for l1 in ref.splitlines():
    if l1.startswith('>'):
        rn1 = l1.strip()[1:]
    else:
        rs1 += stripsequpper(l1.strip())
A1 = [[0.0,0.0,0.0,0.0,0.0] for i in range(len(rs1))]
I1 = [{} for i in range(len(rs1))]

s = ''
dq = 1
ds = 1
for a1 in a.splitlines():
    if 'Plus/Minus' in a1:
        dq=1
        ds=-1
    elif 'Plus/Plus' in a1:
        dq=1
        ds=1
    elif 'Minus/Plus' in a1:
        dq=-1
        ds=1
    elif 'Minus/Minus' in a1:
        dq=-1
        ds=-1
        
    if a1.startswith('Query '):
        QLine1 = a1
    if a1.startswith('Sbjct'):
        SLine1 = a1
        if SwapQuerySubject1:
            QLine1,SLine1 = SLine1,QLine1        
        pq = int(QLine1.split()[1])-1
        q = QLine1.split()[2]

        s = SLine1.split()[2]
        ps = int(SLine1.split()[1])-1
        ins = ''
        for bq,bs in zip(q,s):
            if bq!='-':
                diff1 = 1.0
                diff1 += PosBias1(ps)
                if ps<len(StartBias1):
                    diff1 += StartBias1[ps]
                if ps>=l1-len(EndBias1):
                    diff1 += EndBias1[l1-ps-1]
                if bs.lower()==bq.lower():
                    diff1 += ReferenceBias1
                if nom1>1:
                    diff1 += MultalignBias1
                if ins:
                    if not ins in I1[pq]:
                        I1[pq][ins] = 0.0
                    I1[pq][ins] += diff1
                    ins=''
                A1[pq][BaseD1[bs]] += diff1
                pq += dq
            else:
                ins+=bs
            if bs!='-':
                ps += ds
    if a1.startswith('Range'):
        r1 =int(a1.split()[1].strip(':'))
    if a1.startswith('Sequence ID'):
        n1 = a1.split('Sequence ID:')[1].split('Length')[0] ## name
        l1 = int(a1.split('Length:')[1].split('Num')[0])
        nom1 = int(a1.split()[-1]) ## number of homology segments
c1 = ''
c2 = ''
def Conflicts1(Array4):
    NoConflict1 = 0
    Conflict1 = 0
    if sum(Array4)<1.5:
        return 0,0
    for i in range(4):
        for j in range(4):
            q1 = round(Array4[i])
            q2 = round(Array4[j])
            if i!=j:
                Conflict1+=q1*q2
            if i==j:
                NoConflict1+=(q1*q1-1)
    return NoConflict1,Conflict1

C11, NC11, SI11 = 0,0,0
        
for p0,(a1,i1) in enumerate(zip(A1,I1)):
    vm=0
    if i1:
        vm = max(i1.values())
        lm = sum(i1.values())-vm
    for i11 in sorted(list(i1.keys())):  ##take the first alphabetical insert with valmax instances
        if i1[i11]==vm:
            if lm==0:
                c1 += i11.upper()
            else:
                c1 += i11.lower()
            break
    if sum(a1)==0:
        c1+='n'  ## no data
        c2+='n'  ## no data
    else:
        nc11,c11 = Conflicts1(a1)
        if (nc11,c11)==(0,0):
            SI11+=1.0
        else:
            C11+=c11
            NC11+=nc11
        M = sorted(a1)
        Top = M[-1]
        Next = M[-2]
        AMax = a1.index(Top)
        if Top==Next:
            c1 += 'N'  ## conflicting data -- unresolved [rare]
            c2 += 'N'  ## conflicting data -- unresolved [rare]
        elif Next==0.0:
            c1 += 'GATCD'[AMax]
            c2 += 'GATCD'[AMax]
        else:
            c1 += 'gatcd'[AMax]
            c2 += 'gatcd'[AMax]
        
print '*** Scaffolded sequence with insertions ***'
print c1
print
print '*** Scaffolded sequence with no insertions ***'
print c2
print
print 'Overlap Quality'
print 'SingleRead Positions:',int(SI11)
print 'Disagreeing Read Comparisons: (unit is bases**2)',int(C11)
print 'Agreeing Read Comparisons: (unit is bases**2)',int(NC11)
from time import strftime,localtime
now1 = strftime("%m_%d_%y_%H_%M_%S",localtime())

WithFile1 = 'ScaffoldWithInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+RefFile.split('.')[0].replace('_','')+'_'+now1+'.fasta'
WithoutFile1 = 'ScaffoldWithoutInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+RefFile.split('.')[0].replace('_','')+'_'+now1+'.fasta'
if not(rn1):
    SeqName1 = 'ScaffoldWithInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+RefFile.split('.')[0].replace('_','')
    SeqName2 = 'ScaffoldWithoutInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+RefFile.split('.')[0].replace('_','')
else:
    SeqName1 = 'ScaffoldWithInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+rn1.replace('_','')
    SeqName2 = 'ScaffoldWithoutInsertions_'+AlignFile.split('.')[0].replace('_','')+'_'+rn1.replace('_','')
WF1 = open(WithFile1,mode='w')
WF2 = open(WithoutFile1,mode='w')
WF1.write('>'+SeqName1+LineEnd1)
WF2.write('>'+SeqName2+LineEnd1)
if ReplaceDsWithDashes1:
    c1 = c1.replace('d','-').replace('D','-')
    c2 = c2.replace('d','-').replace('D','-')
for i in range(0,len(c1),LineLengthOut1):
    WF1.write(c1[i:i+100]+LineEnd1)
for i in range(0,len(c2),LineLengthOut1):
    WF2.write(c2[i:i+100]+LineEnd1)
WF1.close()
WF2.close()



