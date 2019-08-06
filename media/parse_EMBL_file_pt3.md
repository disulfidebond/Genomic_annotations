### Description of Workflow

#### Step 3: Create a fasta file of exons that extends the window for mapping reads N base pairs in both directions, where N == the number of base pairs the for the sequencing length.  For Illumina, this is usually 75 or 150.
* Required parameters:
  * ntype -> the type of input, possible options are:
    * rna -> indicates that the input sequenced read data is from rna/cDNA (mutually exclusive with dna)
    * dna -> indicates that the input sequenced read data is from dna (mutually exclusive with rna, do *not* enter dna for cDNA data)
  * filteredfile -> the file containing the input sequenced reads that will be used.  Must be in MHC.dat format
* Other parameters:
  * flength -> the length of padding for the reads, the default is 75
  
                #!/usr/bin/python3
                import argparse
                import time

                def parsedAndCastRange(xList):
                    # requires: string in the format '1..5', returns [1,5] as list of ints
                    unparsedRangeAsList = xList.split('.')
                    parsedRange = list(filter(lambda x: x!= '', unparsedRangeAsList))
                    return [int(parsedRange[0]), int(parsedRange[1])]

                def parsedExonOutput(s,idx):
                    fwdLen = idx
                    revLen = -1*idx
                    retList = []
                    retListExtended = []
                    for i in s:
                        listOfValues = i
                        exonSeq = i[-1]
                        exon_ct = []
                        exon_ct_list = []
                        mxLen = len(listOfValues[3])
                        for j in listOfValues[3]:
                            exon_ct.append(j[0])
                            exon_ct_list.append((j[1][0], j[1][1]))
                        # indexTransition = 0
                        indexTransition = 0
                        checkForIntron = 0
                        eList = []
                        for x in range(0, mxLen, 1):
                            exonStartAndStop = exon_ct_list[x]
                            exonStart = exonStartAndStop[0]
                            if x == 0:
                                exonStart = exonStartAndStop[0] - 1
                            exonStop = exonStartAndStop[1] + 1
                            exonSequence = exonSeq[exonStart:exonStop]
                            indexTransition = exonStartAndStop[1]+1
                            # r = str(listOfValues[0]) + ':' + str(exon_ct[x]) + ',' + exonSequence
                            r = str(listOfValues[0]) + ':' + str(exon_ct[x])
                            eList.append((r, exonSequence))
                        retList.append(eList)
                    for i in range(0, len(retList)):
                        iLength = len(retList[i])
                        seqIntAsString, seqStringList = zip(*retList[i])
                        stringOfAllSeqs = ''.join(seqStringList)
                        rList = []
                        for j in range(0, iLength):
                            currentStringFromList = seqStringList[j]
                            seqID = seqIntAsString[j]
                            if j == 0:
                                cLoc = len(currentStringFromList)
                                endLoc = cLoc + fwdLen
                                endLoc = endLoc - 1
                                paddedString = stringOfAllSeqs[0:endLoc]
                                rList.append((seqID,paddedString))
                            elif j == iLength-1:
                                endString = seqStringList[j-1] + seqStringList[j]
                                paddedString_I = stringOfAllSeqs.find(endString) - fwdLen
                                startPos = paddedString_I - fwdLen
                                paddedStringIdx = stringOfAllSeqs.find(endString) + len(endString)
                                paddedStringIdx = paddedStringIdx
                                paddedString = stringOfAllSeqs[paddedString_I:paddedStringIdx]
                                rList.append((seqID,paddedString))
                            else:
                                paddedString_Ipos = stringOfAllSeqs.find(str(seqStringList[j-1]))
                                startCheck = paddedString_Ipos - fwdLen
                                fwdLenPadded = stringOfAllSeqs.find(str(seqStringList[j+1])) + fwdLen
                                if len(seqStringList[j+1]) < 10:
                                    endString = seqStringList[j-1] + seqStringList[j]
                                    fwdLenPadded = stringOfAllSeqs.find(endString) + len(endString) + fwdLen
                                # stopCheck = fwdLenPadded - len(stringOfAllSeqs)
                                startPos = startCheck
                                stopPos = fwdLenPadded
                                # scenario 1: fwdLen > strings before currentString
                                if startCheck < 0:
                                    startPos = 0
                                # scenario 2: fwdLen > strings after currentString
                                if fwdLenPadded > len(stringOfAllSeqs):
                                    stopPos = len(stringOfAllSeqs)
                                    # NOTE: this is a bit of a hack, and assumes all sequences after index 3 will have decreasing length
                                    # create a better one when not as sleepy
                                # scenario 3 pos(prevString) < fwdLen < pos(nextString): keep values unchanged
                                paddedString = stringOfAllSeqs[startPos:stopPos]
                                rList.append((seqID,paddedString))
                        retListExtended.append(rList)
                    return retListExtended

                if __name__ == "__main__":
                    parser = argparse.ArgumentParser(description='parser for MHC file from EMBL')
                    parser.add_argument('-f', '--filteredfile', help='filtered file from EMBL', required=True)
                    parser.add_argument('ntype', choices=['rna', 'dna'], help='Enter either \'rna\' for rna/cDNA input, or \'dna\' for dna input')
                    parser.add_argument('-i', '--flength', help='read length', required=True)
                    args = parser.parse_args()


                    alleleList = []
                    with open(args.filteredfile, 'r') as f:
                        s_id = ''
                        s_alleleID = ''
                        s_exonList = []
                        alleleExonRange = []
                        alleleExonSeq = []
                        exonSeqs = ''
                        exonCt = 1
                        seqFlag = False
                        for i in f:
                            i = i.rstrip('\r\n')
                            sLine = i.split(' ')
                            sLineList = list(filter(lambda x: x != '', sLine))
                            if i[0:2] == 'KW':
                                s_id = sLineList[1]
                            elif i[0:2] == 'ID':
                                s_alleleID = sLineList[1]
                            elif i[0:2] == 'FT':
                                if sLineList[1] == 'allele':
                                    alleleExonRange = parsedAndCastRange(sLineList[2])
                                elif sLineList[1] == 'exon':
                                    exonSeqRange = parsedAndCastRange(sLineList[2])
                                    exonN = 'exon' + str(exonCt)
                                    alleleExonSeq.append((exonN, exonSeqRange))
                                    exonCt += 1
                                else:
                                    continue
                            elif i[0:2] == 'SQ':
                                seqFlag = True
                                continue
                            else:
                                if seqFlag:
                                    if i[0:2] == '//':
                                        seqFlag = False
                                        alleleList.append([s_id, s_alleleID, alleleExonRange, alleleExonSeq, exonSeqs])
                                        alleleExonSeq = []
                                        exonCt = 1
                                        exonSeqs = ''
                                        continue
                                    else:
                                        s = ''.join(sLineList[:-1])
                                        exonSeqs += s
                                else:
                                    continue
                                    
                    if args.outputType == 'exonList':
                        res = parsedExonOutput(alleleList, args.flength)
                        for i in res:
                            print(i)
                    else:
                        res = parsedExonOutput(alleleList, args.flength)
                        for i in res:
                            for j in i:
                                n = i[0]
                                v = i[1]
                                hString = '>' + str(n)
                                seqString = str(v)
                                print(hString + '\n' + seqString)
