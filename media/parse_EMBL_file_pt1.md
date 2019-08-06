### Description of Workflow

#### Step 1. download MHC dataset from IPD.  It will have the file extension ".dat", however it can be viewed/edited as a text file.  Note that it is very large, so opening it in Atom, BBEdit, or similar text editors/IDE is strongly discouraged.

A preview of this dataset is available here:

                    login$ head -n 100 MHC_dat.txt 
                    ID   NHP00001
                    XX   
                    DT   15/07/2008 (Release)
                    XX   
                    KW   Aona-DQA1*27:01
                    XX   
                    DR   EMBL; AF201293.
                    XX   
                    CC   The nucleotide sequence provided is a CDS sequence, constructed from the 
                    CC   sequences submitted to the IPD-MHC Database. The sequence below is the 
                    CC   official sequence for Aona-DQA1*27:01 as approved by the MHC Nomenclature 
                    CC   Committee and as a result the sequence described in any cross references 
                    CC   may differ from that shown in the IPD-MHC 
                    CC   Database.
                    XX   
                    OS   Aotus nancymaae (Nancy Ma's Night Monkey)
                    XX   
                    OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; 
                    OC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Platyrrhini; 
                    OC   Aotidae; Aotus; 
                    XX   
                    FH   Key            Location/Qualifier
                    FH   
                    FT   allele         1..226
                    FT                  /status="public"
                    FT   source         1..226
                    FT                  /db_xref="taxon:37293"
                    FT   exon           1..226
                    FT                  /number=2
                    FT   CDS            1..226
                    FT                  /gene="DQA1"
                    FT                  /allele="Aona-DQA1*27:01"
                    FT                  /translation="DHVAAYGINLYQSYGLSGQYTHEFDGDEEFYVDLGRKETVWRLPVF
                    FT                  SKFAGFDPQGALTNIAAGKHNLDILIKR"
                    FT                  /codon_start=3
                    XX   
                    SQ   Sequence 226 BP; 57 A; 53 C; 58 G; 58 T; 0 other;
                         CTGACCATGT TGCCGCTTAC GGTATAAACT TGTACCAGTC TTATGGTCTC TCTGGCCAGT        60
                         ACACCCACGA ATTTGATGGA GATGAGGAGT TCTACGTGGA CCTGGGAAGA AAGGAGACTG       120
                         TCTGGCGATT GCCTGTGTTC AGCAAATTTG CAGGTTTTGA CCCTCAGGGT GCACTGACAA       180
                         ACATCGCTGC GGGAAAACAC AACTTGGACA TCCTGATTAA ACGCTC                      226
                    //
                    ID   NHP00002
                    XX   
                    DT   15/07/2008 (Release)
                    XX   
                    KW   Aona-DQA1*27:02
                    XX   
                    DR   EMBL; AF201294.
                    XX   
                    CC   The nucleotide sequence provided is a CDS sequence, constructed from the 
                    CC   sequences submitted to the IPD-MHC Database. The sequence below is the 
                    CC   official sequence for Aona-DQA1*27:02 as approved by the MHC Nomenclature 
                    CC   Committee and as a result the sequence described in any cross references 
                    CC   may differ from that shown in the IPD-MHC 
                    CC   Database.
                    XX   
                    OS   Aotus nancymaae (Nancy Ma's Night Monkey)
                    XX   
                    OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; 
                    OC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Platyrrhini; 
                    OC   Aotidae; Aotus; 
                    XX   
                    FH   Key            Location/Qualifier
                    FH   
                    FT   allele         1..226
                    FT                  /status="public"
                    FT   source         1..226
                    FT                  /db_xref="taxon:37293"
                    FT   exon           1..226
                    FT                  /number=2
                    FT   CDS            1..226
                    FT                  /gene="DQA1"
                    FT                  /allele="Aona-DQA1*27:02"
                    FT                  /translation="DHVAAYGINLYQSYGLSGQYTHEFDGDEEFYMDLERKETVWRLPVF
                    FT                  SKFAGFDPQGALTNIAAGKHNLDILIKR"
                    FT                  /codon_start=3
                    XX   
                    SQ   Sequence 226 BP; 59 A; 53 C; 56 G; 58 T; 0 other;
                         CTGACCATGT TGCCGCTTAC GGTATAAACT TGTACCAGTC TTATGGTCTC TCTGGCCAGT        60
                         ACACCCACGA ATTTGATGGA GATGAGGAGT TCTACATGGA CCTGGAAAGA AAGGAGACTG       120
                         TCTGGCGATT GCCTGTGTTC AGCAAATTTG CAGGTTTTGA CCCTCAGGGT GCACTGACAA       180
                         ACATCGCTGC GGGAAAACAC AACTTGGACA TCCTGATTAA ACGCTC                      226
                    //

#### Step 2: Parse out the desired values.  Several methods are available, the example shown here is usable but strongly in need of revision and automation.  The output file, here named 'parsed_mhc_output.txt', will have all of the entries for macaca mulatta (or whatever)

                # parse out all entries with the string 'Macaca mulatta', this will create a noisy list
                grep -B16 'Macaca fascicularis' MHC.dat > macaca_mulatta_listOfIDs.mhc.txt
                # these two commands remove the noise and provides a list of ID's only
                grep 'ID   ' macaca_fascicularis_listOfIDs.mhc.txt > macaca_fascicularis_listOfIDs.mhc.parsed.txt 
                perl -p -e 's/ID\s+//g' macaca_fascicularis_listOfIDs.mhc.parsed.txt > macaca_fascicularis_listOfIDs.mhc.parsedID.txt
                
                # then run this short python script, output is to STDOUT
                
                #!/usr/bin/python
                import time

                def importEntries(f):
                    pList = []
                    with open(f) as fOpen:
                        currentItem = []
                        for i in fOpen:
                            i = i.rstrip('\r\n')
                            if i[0:3] == '//':
                                pList.append(currentItem)
                                currentItem = []
                            else:
                                currentItem.append(i)
                    return pList

                def parseOutSelectedEntries(l, itmList):
                    rList = []
                    for itm in l:
                        checkItem = itm[0]
                        splitItem = checkItem.split(' ')
                        checkItemFilteredList = list(filter(lambda x: x != '', splitItem))
                        if checkItemFilteredList[1] in itmList:
                            rList.append(itm)
                    return rList

                idList = []
                with open('macaca_fascicularis_listOfIDs.mhc.parsedID.txt') as f:
                    for i in f:
                        i = i.rstrip('\n\r')
                        idList.append(i)
                importedList = importEntries('MHC.dat.txt')
                parsedAndImportedList = parseOutSelectedEntries(importedList, idList)
                for itm in parsedAndImportedList:
                    o = '\n'.join(itm)
                    o += '\n##\n'
                    # some flag is needed for the next parser step.
                    # An alternative is to keep the '//' above and use that instead
                    # For reasons that are not entirely clear, do *NOT* use a string longer than 2 characters
                    print(o)

