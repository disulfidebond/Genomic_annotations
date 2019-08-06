### Description of Workflow

#### Step 2: Use python file to parse the filtered output from part 1 (most likely the file named 'parsed_mhc_output.txt'). The python script uses argparse(), so it can be called from the commandline with one of the following options, but note that output is to STDOUT, so it will need to be redirected to a file.
* exonList -> a text file with a comma separated list of exons
* fastaExonList -> a fasta formatted list of exons
* mergedFastaList -> a fasta formatted list of cDNA sequences

                #!/usr/bin/python3
                import argparse
                import time

                def parsedAndCastRange(xList):
                    # requires: string in the format '1..5', returns [1,5] as list of ints
                    unparsedRangeAsList = xList.split('.')
                    parsedRange = list(filter(lambda x: x!= '', unparsedRangeAsList))
                    return [int(parsedRange[0]), int(parsedRange[1])]

                def parsedExonOutput(s):
                    retList = []
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
                        for x in range(0, mxLen, 1):
                            exonStartAndStop = exon_ct_list[x]
                            exonStart = exonStartAndStop[0]
                            if x == 0:
                                exonStart = exonStartAndStop[0] - 1
                            exonStop = exonStartAndStop[1] + 1
                            exonSequence = exonSeq[exonStart:exonStop]
                            indexTransition = exonStartAndStop[1]+1
                            r = str(listOfValues[0]) + ':' + str(exon_ct[x]) + ',' + exonSequence
                            retList.append(r)
                    return retList


                if __name__ == "__main__":
                    parser = argparse.ArgumentParser(description='parser for MHC file from EMBL')
                    parser.add_argument('-f', '--filteredfile', help='filtered file from EMBL', required=True)
                    parser.add_argument('outputType', choices=['exonList', 'fastaExonList', 'mergedFastaList'], help='Type one of three output options without quotes: \'exonList\' outputs a list with Allele:exon,exonsequence ; \'fastaExonList\' outputs a fasta file split into exons ; \'mergedFastaList\' outputs a fasta file with the merged exons for each allele and the header name of the allele')
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
                            if i[0:3] == 'KW':
                                s_id = sLineList[1]
                            elif i[0:3] == 'ID':
                                s_alleleID = sLineList[1]
                            elif i[0:3] == 'FT':
                                if sLineList[1] == 'allele':
                                    alleleExonRange = parsedAndCastRange(sLineList[2])
                                elif sLineList[1] == 'exon':
                                    exonSeqRange = parsedAndCastRange(sLineList[2])
                                    exonN = 'exon' + str(exonCt)
                                    alleleExonSeq.append((exonN, exonSeqRange))
                                    exonCt += 1
                                else:
                                    continue
                            elif i[0:3] == 'SQ':
                                seqFlag = True
                                continue
                            else:
                                if seqFlag:
                                    if i[0:2] == '##':
                                        # modify this as necessary.  It is unknown why python string indices 
                                        # are completely wacky--possibly something Bash does?
                                        # be advised you may need to change the range from 0:3 
                                        # or change the string range in the parsing above to 0:2 
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
                        res = parsedExonOutput(alleleList)
                        for i in res:
                            print(i)
                    else:
                        if args.outputType == 'mergedFastaList':
                            for i in alleleList:
                                hString = '>' + str(i[0]) + ' ' + str(i[1])
                                seqString = str(i[-1])
                                print(hString + '\n' + seqString)
                        else:
                            res = parsedExonOutput(alleleList)
                            for i in res:
                                iSplit = i.split(',')
                                hString = '>' + str(iSplit[0])
                                seqString = str(iSplit[1])
                                print(hString + '\n' + seqString)
