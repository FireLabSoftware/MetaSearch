import matplotlib.pyplot as plt
import numpy as np
import cPickle
import sys

def revComp(seq):
    seq=seq.upper()
    seq=seq.replace('A','t')
    seq=seq.replace('T','a')
    seq=seq.replace('C','g')
    seq=seq.replace('G','c')
    seq=seq.upper()
    seq=seq[::-1]
    return seq

def dictIncrement(seq,dictionary):
    if seq in dictionary:
        dictionary[seq]+=1
    else:
        dictionary[seq]=1

def dictSorter(dictionary):
    return sorted(dictionary,key=lambda x: dictionary[x], reverse=True)

def checkNoN(seq):
    seqToCheck=seq.upper()
    if 'N' in seqToCheck:
        return False
    else:
        return True

class alignedObject:
    def __init__(self,alignmentInfo,refUsed,count):
        self.alignmentInfo=alignmentInfo
        self.refUsed=refUsed
        self.count=count

class alignmentInfoEvent:
    def __init__(self,typeOfEvent,position,notes):
        self.typeOfEvent=typeOfEvent
        self.position=position
        self.notes=notes

def printAlignmentInfoList(alignmentInfoList):
    listToPrint=[]
    for eachEntry in alignmentInfoList:
        listToPrint.append((eachEntry.typeOfEvent,eachEntry.position,eachEntry.notes))
    print listToPrint       

#when n=0, if a sequence starts with a diversity base, then behaviour of alignment was to treat that base as an insertion at beginning
#thus made n=1 below

def smithWaterman(seq1,seq2,m=2,s=-1,d=-2,n=1,collapseInsertions=True): 
    #m is for matches, s for mismatches, d for gaps, n is for ns in reference--sites where diversity there
    seq1=seq1.upper()
    seq2=seq2.upper()
    dpTable=np.zeros((len(seq1)+1,len(seq2)+1)) #assume seq1 is reference
    pointerTable=[]
    seqEmitTable=[]
    dummyVector=[]
    for i in range(0,len(seq2)+1):
        dummyVector.append(0)
    for j in range(0,len(seq1)+1):
        pointerTable.append(list(dummyVector))
        seqEmitTable.append(list(dummyVector))
    #initialization--not required as everything already set to 0, this is a big change from needlemanWunsch
    #update dynamic programming table
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if seq1[i-1]=='N':
                subScore=n
            else:
                if seq1[i-1]==seq2[j-1]:
                    subScore=m
                else:
                    subScore=s
            diagScore=dpTable[i-1][j-1]+subScore
            topDownScore=dpTable[i-1][j]+d
            leftRightScore=dpTable[i][j-1]+d
            #add constraint compared to needlemanWunsch that minimum value of each cell is 0 to get smith-waterman
            if diagScore>=topDownScore and diagScore>=leftRightScore and diagScore>0:
                dpTable[i][j]=diagScore
                pointerTable[i][j]=(i-1,j-1)
                seqEmitTable[i][j]=(seq1[i-1],seq2[j-1],i-1,j-1)
            elif topDownScore>=diagScore and topDownScore>=leftRightScore and topDownScore>0:
                dpTable[i][j]=topDownScore
                pointerTable[i][j]=(i-1,j)
                seqEmitTable[i][j]=(seq1[i-1],'_',i-1,'deletion')
            elif leftRightScore>=diagScore and leftRightScore>=topDownScore and leftRightScore>0:
                dpTable[i][j]=leftRightScore
                pointerTable[i][j]=(i,j-1)
                seqEmitTable[i][j]=('_',seq2[j-1],'insertion',j-1)
            elif diagScore<=0 and topDownScore<=0 and leftRightScore<=0:
                dpTable[i][j]=0
                pointerTable[i][j]=0
                seqEmitTable[i][j]=0

    #construct alignment--first find max value in dpTable
    maxIndexTuple=np.unravel_index(np.argmax(dpTable),dpTable.shape)                # using convention of finding maximum to resolve any ties in alignment, can refine this later to deal with non-unique alignments, also perhaps later want to see the score for the second best alignment
    #print dpTable
    #print pointerTable
    '''
    if dpTable[maxIndexTuple[0]][maxIndexTuple[1]]<=0:
        #have a condition to test for this--not required as at least one base should align and hence score should be >0
        return
    '''
    #dpTable[maxIndexTuple[0]][maxIndexTuple[1]] should be alignment score
    #alignmentList=[]
    alignmentInfo=[]
    currentCoordinate=maxIndexTuple
    lastPoint=1
    while dpTable[currentCoordinate[0]][currentCoordinate[1]]!=0:
        relevantSeqEmit=seqEmitTable[currentCoordinate[0]][currentCoordinate[1]]
        #alignmentList.append(relevantSeqEmit)
        if lastPoint==1:
            alignmentInfo.append(alignmentInfoEvent('End',currentCoordinate[0],'')) #end cannot be at an insertion or deletion as it costs points
            if currentCoordinate[1]==len(seq2): #no insertion there at 3' end
                pass
            else:
                alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],seq2[currentCoordinate[1]:]))            
            lastPoint=0
        if relevantSeqEmit[3]=='deletion':
            alignmentInfo.append(alignmentInfoEvent('Deletion',currentCoordinate[0],''))
        elif relevantSeqEmit[2]=='insertion':
            alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],relevantSeqEmit[1]))
        else:
            if relevantSeqEmit[0]=='N':
                alignmentInfo.append(alignmentInfoEvent('Diversity',currentCoordinate[0],relevantSeqEmit[1]))
            elif relevantSeqEmit[0]==relevantSeqEmit[1]:
                alignmentInfo.append(alignmentInfoEvent('Match',currentCoordinate[0],relevantSeqEmit[1]))
            else:
                alignmentInfo.append(alignmentInfoEvent('Mismatch',currentCoordinate[0],relevantSeqEmit[1]))
        currentCoordinate=pointerTable[currentCoordinate[0]][currentCoordinate[1]]
    if currentCoordinate[1]==0: #no insertion there at 5' end
        pass
    else:
        alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],seq2[:currentCoordinate[1]]))
    alignmentInfo.append(alignmentInfoEvent('Start',currentCoordinate[0],'')) 

    alignmentInfo=alignmentInfo[::-1]

    if collapseInsertions==True:
        #collapse insertions that occur together into one sequence    
        newAlignmentInfo=[]
        insertionFound=0
        insertCoord=0
        insertSeq=''    
        for eachEl in alignmentInfo:
            if eachEl.typeOfEvent=='Insertion':            
                insertionFound=1
                insertCoord=eachEl.position
                insertSeq=insertSeq+eachEl.notes
            else:
                if insertionFound==1:
                    newAlignmentInfo.append(alignmentInfoEvent('Insertion',insertCoord,insertSeq))
                    insertCoord=0
                    insertSeq=''
                    insertionFound=0
                newAlignmentInfo.append(eachEl)                    
        return [dpTable[maxIndexTuple[0]][maxIndexTuple[1]],newAlignmentInfo]
    else:
        return [dpTable[maxIndexTuple[0]][maxIndexTuple[1]],alignmentInfo]

