## Overview
The following functions take as input:
* alleleString -> string of nucleotides that will attempt to be matched to a known exon string
* alleleName -> string that is the name linked to alleleString, which will first be used to match to a known exon string
* alleleExonTuple -> Exons must be split into a list of tuples with the format:

        [
          ('Allele1Name', ['exon1-string', 'exon2-string',...,'exonN-string']),
          ('Allele2Name', ['exon1-string', 'exon2-string',...,'exonN-string']),
          ('Allele3Name', ['exon1-string', 'exon2-string',...,'exonN-string'])
        ]
        # example:
        ('Mamu-A1*001:01:01:01', ['CCCC', 'CCCCCCCC', 'CGGAA'])

* alleleRgx -> regular expression pattern that will be used to split the alleleName, typically ':'

The functions work as follows:

1) Attempt to match alleleName to an existing name. If an exact match succeeds, proceed to step 3. Otherwise, continue to step 2.

2) Split allleName and walk the naming tree.  For example, the allele Mamu-A1*001:01:01:01 would first attempt to match any existing alleles named Mamu-A1*001:01:01, then Mamu-A1*001:01, then Mamu-A1*001. If no match is found, throw error and exit. Otherwise, proceed to step 3.

3) For each allele in alleleExonTuple, attempt to match each exon in order to the alleleString. If a non-overlapping match is found, return the positions of the match as (startPos, endPos). If a match is found that is overlapping, return the start position and -1 as (startPos, -1). Otherwise, return (-1,-1).

4) Check the results. If there is more than one instance of an overlap or a match failure, set an integer bool flag. Return the result as:

`(0 || 1, alleleName, matched-AlleleName, [alleleString, [matched-AlleleExonStrings], [alleleString-matchedPos]])`

## Code

        def exonPosMatching(startPosStart, seq, exonString):
            startPos = seq.find(exonString)
            if startPos < startPosStart:
                return (startPos, -1)
            elif startPos > -1:
                endPos = startPos + len(exonString)
                return (startPos, endPos)
            else:
                return (-1, -1)

        def matchExonsToString(alleleString, alleleName, alleleExonTuple, alleleRgx):
            res_pos = []
            _start_pos = 0
            filter_res = list(filter(lambda x: x[0] == alleleName, alleleExonTuple))
            if filter_res:
                _filter_res = filter_res[0][1][:-1]
                filter_res_lastExon = filter_res[0][1][-1]
                for i in _filter_res[0][1]:
                    matched = exonPosMatching(_start_pos, alleleString, i)
                    _start_pos = matched[0]
                    res_pos.append(matched)
            else:
                string_scan = alleleName.split(alleleRgx)
                stopCt = len(string_scan) # + 1
                for i in range(1, stopCt):
                    res = alleleRgx.join(string_scan[:-i])
                    filter_res = list(filter(lambda x: x[0] == res, alleleExonTuple))
                    if filter_res:
                        _filter_res = filter_res[0][1][:-1]
                        filter_res_lastExon = filter_res[0][1][-1]
                        for x in _filter_res:
                            matched = exonPosMatching(_start_pos, s_seq, x)
                            _start_pos = matched[0]
                            res_pos.append(matched)
                        break # will stop after first match
            if not res_pos:
                return None
            end_pos = _start_pos + len(filter_res[0][1][-2])
            lastExonPosSearch = s_seq[end_pos:]
            matched = exonPosMatching(0, lastExonPosSearch, filter_res_lastExon)
            if matched[1] > -1:
                matched1 = matched[0] + end_pos
                matched2 = matched[1] + end_pos
                matched = (matched1, matched2)
            res_pos.append(matched)
            checkPos = list(filter(lambda x: x[1] == -1, res_pos))
            if len(checkPos) > 1:
                return (0, alleleName, filter_res[0][0], [alleleString, filter_res[0][1], res_pos])
            else:
                return (1, alleleName, filter_res[0][0], [alleleString, filter_res[0][1], res_pos])
