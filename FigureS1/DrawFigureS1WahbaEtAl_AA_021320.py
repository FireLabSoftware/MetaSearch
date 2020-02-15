## This is a dedicated script to draw Figure S1 in Wahba et al, as submitted 02.13.20
## The script passes a set of objects for drawing to the VSG module for realization on a 2D canvas
## Details of the analysis are from Wahba et al, 2020 (submitted).  VSG_ModuleDB is available at github/FireLabSoftware
## Copywrite 2020 Fire Lab and Stanford University, with no guarantees

from sys import argv
from VSG_ModuleDB import *
from operator import ne
try:
    from itertools import imap
except:
    imap=map


RefFile = 'NC_045512.2.fasta'
Starts1 = (266, 13468, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558)
Ends1 = (13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674)
ProteinList1 = ((266,805,"leaderP"),
        (806,2719,"nsp2"),
        (2720,8554,"nsp3"),
        (8555,10054,"nsp4"),
        (10055,10972,"3C-like-protease"),
        (10973,11842,"nsp6"),
        (11843,12091,"nsp7"),
        (12092,12685,"nsp8"),
        (12686,13024,"nsp9"),
        (13025,13441,"nsp10"),
        (13442,13480,"nsp11"),
        (13442,16236,"RdRP"),
        (16237,18039,"helicase"),
        (18040,19620,"3'-to-5'exonuclease"),
        (19621,20658,"endoRNAse"),
        (20659,21552,"2'O-ribose MeTransferase"),
        (21563,25384,"Spike"),
        (25393,26220,"ORF3a"),
        (26245,26472,"E"),
        (26523,27191,"M"),
        (27202,27387,"ORF6"),
        (27394,27759,"ORF7a"),
        (27756,27887,"ORF7b"),
        (27894,28259,"ORF8"),
        (28274,29533,"N"),
        (29558,29674,"ORF10"))

ScaffoldFileD1 ={'ScaffoldWithoutInsertions_41X3HN2V01R-Alignment_NC045512_02_10_20_07_44_58.fasta':('magenta','red','orange'),
                'ScaffoldWithoutInsertions_BlastALignRaTG13Wuhan_NC045512_02_11_20_08_00_22.fasta':('cyan','blue','green')}
ScaffoldFileD2 ={'ScaffoldWithoutInsertions_41X3HN2V01R-Alignment_NC045512_02_10_20_07_44_58.fasta':('SRR10168377 Read-Aggregate [Pangolin]'),
                'ScaffoldWithoutInsertions_BlastALignRaTG13Wuhan_NC045512_02_11_20_08_00_22.fasta':('RaTG13 [Bat]')}
## Input files
if len(argv)>=3:
    RefFile = argv[1]  ## Original Reference File for the search in FastA format (one sequence only please)
    ScaffoldFile = argv[2] ## File of alignments from NCBI blast output on web

Window1 = 75
hScale1 = 5.0
vScale1 = 150.0
cRad1 = 3.0
LineWidth1 = 1.5


## a fast genetic code parser.  Everything needs to be upper case and needs to be a real base
## can pre-filter also, or at least make sure of upper case with translate(s.upper()), etc
a1list=[]
for a1 in ['A','G','C','T']:
    for a2 in ['A','G','C','T']:
        for a3 in ['A','G','C','T']:
            a1list.append(a1+a2+a3)
aa='KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF'
inferExtend={'AC':'T','GG':'G','GC':'A','GT':'V',
             'CG':'R','CC':'P','CT':'L','TC':'S'}
gc1={a1list[inde]:aa[inde] for inde in range(64)}
def Tr1(s):
    if s in gc1:
        return gc1[s]
    else:
        return ''

def translate(s):
    lt=(len(s)//3)*3   
    if s[lt:] in inferExtend:
        sg=s+'G'  ## add an arbitrary base, since for these 2/3 codons the amino acid is dertermines
        return ''.join(gc1[sg[i:i+3]] for i in range(0,lt+3,3))
    else:
        return ''.join(gc1[s[i:i+3]] for i in range(0,lt,3))

def ham1(a,b):
    return(sum(imap(ne,a,b))+abs(len(a)-len(b)))  ## quick hamming distance (with penalties for any mismatch in length


vset(bg=white)
ref = open(RefFile,mode='rU').read()
rs1 = ''
rn1 = ''
for l1 in ref.splitlines():
    if l1.startswith('>'):
        rn1 = l1.strip()[1:]
    else:
        rs1 += l1.strip()
ls1 = len(rs1)

P1 = [-1]*ls1
Q1 = [-1]*ls1
for ts1,es1 in zip(Starts1,Ends1):
    for tp1 in range(ts1-1,es1):
        P1[tp1]=(tp1-ts1+1)%3
        Q1[tp1]=tp1-P1[tp1]-1
        
    

for ScaffoldFile1 in ScaffoldFileD1.keys():
    var = open(ScaffoldFile1,mode='rU').read()
    ChartColor1 = ScaffoldFileD1[ScaffoldFile1]
    Source1 = ScaffoldFileD2[ScaffoldFile1]
    vs1 = ''
    vn1 = ''
    for l1 in var.splitlines():
        if l1.startswith('>'):
            vn1 = l1.strip()[1:]
        else:
            vs1 += l1.strip()
    dS1 = [0 for i in range(ls1)]  ## Synonymous differences
    dN1 = [0 for i in range(ls1)]  ## NonSynonymous differences
    dO1 = [0 for i in range(ls1)]  ## Other differences
    for i in range(len(rs1)):
        if vs1[i]!=rs1[i]:
            if Q1[i]==-1 and vs1[i].upper()!='N' and rs1[i].upper()!='N':
                dO1[i]=1
            else:
                Codon1 = Q1[i]
                OldAA1 = Tr1(rs1[Codon1:Codon1+3])
                NewAA1 = Tr1(vs1[Codon1:Codon1+3])
                if OldAA1 or NewAA1:
                    if OldAA1 != NewAA1:
                        dN1[i] = 1
                    else:
                        dS1[i] = 1
    tS1 = [0 for i in range(ls1)]  ## Total Synonymous differences
    tN1 = [0 for i in range(ls1)]  ## Total NonSynonymous differences
    tO1 = [0 for i in range(ls1)]  ## Total Other differences
    tC1 = [0 for i in range(ls1)] ## Total in region
    for i in range(ls1-Window1):
        Chunk1 = vs1[i:i+Window1].upper()
        Chunk2 = rs1[i:i+Window1].upper()
        if not('N' in Chunk1):        
            xS1 = sum(dS1[i:i+Window1])
            xN1 = sum(dN1[i:i+Window1])
            xO1 = sum(dO1[i:i+Window1])
            for j in range(i,i+Window1):
                if 0 <= j < ls1:
                    mult1 = 1.0/(1.0+abs(j-(i+Window1/2.0)))
                    tS1[j]+=xS1*mult1
                    tN1[j]+=xN1*mult1
                    tO1[j]+=xO1*mult1
                    tC1[j]+=1*mult1
    oldx1 = -1             
    for i in range(0,ls1):
        if tC1[i]>0 and vs1[i].upper()!='N' and rs1[i].upper()!='N':
            vvS1 = tS1[i]*1.0/tC1[i]        
            vvN1 = tN1[i]*1.0/tC1[i]        
            vvO1 = tO1[i]*1.0/tC1[i]        
            vcircle(xc=i/hScale1, yc=vvS1*vScale1, r=cRad1, fill=ChartColor1[0], stroke=none, strokwidth=0, xg=i, yg=(100.0*vvS1)/Window1, colorkey='Synonymous: '+Source1)
            vcircle(xc=i/hScale1, yc=vvN1*vScale1, r=cRad1, fill=ChartColor1[1], stroke=none, strokwidth=0, xg=i, yg=(100.0*vvN1)/Window1, colorkey='Non-Synonymous: '+Source1)
            vcircle(xc=i/hScale1, yc=vvO1*vScale1, r=cRad1, fill=ChartColor1[2], stroke=none, strokwidth=0, xg=i, yg=(100.0*vvO1)/Window1, colorkey='Non-Coding: '+Source1)
            if oldx1>=0:
                vline(x1=oldx1,x2=i/hScale1,y1=oldS1,y2=vvS1*vScale1, stroke=ChartColor1[0], strokwidth=LineWidth1)
                vline(x1=oldx1,x2=i/hScale1,y1=oldN1,y2=vvN1*vScale1, stroke=ChartColor1[1], strokwidth=LineWidth1)
                vline(x1=oldx1,x2=i/hScale1,y1=oldO1,y2=vvO1*vScale1, stroke=ChartColor1[2], strokwidth=LineWidth1)
            oldx1 = i/hScale1
            oldS1 = vvS1*vScale1
            oldN1 = vvN1*vScale1
            oldO1 = vvO1*vScale1
        else:
            oldx1 = -1
vgrid(gxlabel="Position in NC045512 [2019-nCoV]", gylabel="Variation Index (%)", gtitle='Inferred Variation: Windowed Average')
yq = -800
for k1,(st1, en1, nm1) in enumerate(ProteinList1):
    vrect(x1=st1/hScale1,x2=en1/hScale1,yc=yq-(k1%8)*80,height=70,stroke=none,strokewidth=0,fill=(127,127,127))
    vtext(text=nm1,x1 = en1/hScale1+10, yc=yq-(k1%8)*80,font="DejaVuBold 96",color=black)
vcolorkey(y1=-2000)
vdisplay("PositionalPlot.svg")


    
    





