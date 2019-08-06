# Genomic annotations and file formats
A recent project involved manipulating and reformating common Bioinformatics file formats, such as Genbank, EMBL, and fasta. The required tasks were only vaguely connected, so this repository is a cookbook for file manipulation of common Bioinformatics file formats, including parsing and I/O. The general outline for each section is to explain the file format, explain possible tasks, and then detail how BioPython can be used to accomplish these goals, and at least one possibility for completing the task(s) from scratch.

## FASTA
The [Fasta](https://zhanglab.ccmb.med.umich.edu/FASTA/) file format is a text file with the following requirements:

        >headerLine
        Sequence-Line
        
1) The headerLine must always begin with '>', and is terminated by a newline.  Spaces and special characters are allowed, but may not be advisable.
2) The Sequence-Line above can **only** have the characters 'A,T,C,G'. Whitespace characters and special characters are never allowed, and a newline character immediately signals the end of the sequence.
  * exception: [IUPAC Nucleotide characters](https://www.bioinformatics.org/sms/iupac.html) are permitted, however, these characters are not accepted by all software programs, and may cause errors without warning or with no explanation.
3) Newline characters are **only** allowed after the header line and the end of the sequence line.
4) Multiple fasta entries are allowed in the same file, provided rules 1,2,3 above are followed.  For example, this is a valid fasta file:

        >headerLine 1
        AATTCCGGAATTCCAA
        >header Line 2
        CCCCCCCAAAAATTTTT
        >HEADERLINE|3
        CCGGTTAAAGGAAAAA
        
Beyond there, standardization ends, and problems usually begin. The following fasta file is completely valid, albeit also complete gibberish

        >FDSAH | I forgot what the header was
        AACCTTGGACACAGATATAGAGCCAACACGGTTAACCAA

### Recommendations
* Use meaningful names in the fasta header line, similar to standard coding conventions
* Use a descriptive file name, and be aware that if the filename does not end in 'fasta', 'fa', or 'fas', it may be rejected by a Bioinformatics software program.

### Code Examples
Creating a fasta file is fairly straightforward. The examples for parsing a fasta file can be modified to create a fasta file from input text. 

To parse a fasta file, you could use something similar to:

        # bash
        grep '>' someFastaFile # outputs the fasta file header(s)

Or via BioPython:

        from Bio import SeqIO
        for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
          print(seq_record.id)
          print(repr(seq_record.seq))
          print(len(seq_record))
        # output would be the header, the sequence characters, and the length

## EMBL
The EMBL file format is very similar to the GenBank file format. In addition to the nucleotide sequence and features, it contains additional information like the identity of the submitter, and publications associated with the nucleotide sequence. It has the general format of a two or threee character identifier for each newline, followed by 3 whitespace characters, followed by additional data. Tab characters are not allowed. [An example EMBL file is included with this repository](https://github.com/disulfidebond/Genomic_annotations/blob/master/media/exampleEMBL.txt), and the following figure describes how the file must be formatted:
![](https://github.com/disulfidebond/Genomic_annotations/blob/master/media/EMBL_image1.png)
A non-exhaustive list of the character identifiers is:

* ID -> The identification line. It **must** have the format:

        # example:
        CAA39891; SV 1; linear; genomic DNA; STD; PRO; 225 BP.
        # format:
        Identifier; sequence_version{N} ; topology{'linear','circular'}; molecule_type={'genomic DNA'}; STD; data_file_division=PRO; total-length-in-BP
        
* AC -> Accession ID
* PR -> Project, typically the [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) identifier
* FH -> Feature Header, this described the features that are listed in the FT line(s)
* FT -> Feature Type, a given feature type, description, or information
* XX -> blank line
* // -> Signals the end of an entry. This must be present, even if the EMBL file has only one entry

### Recommendations
* Pay attention to details. 
* You may have to manually curate some of the entries to ensure that the EMBL file meets the required specifications for the software program that you are using.

### Code Examples
You can use BioPython to parse or create an EMBL file. Note that BioPython will not provide helpful error messages if the EMBL file has syntactical errors.

        for record in SeqIO.parse("SomeFile.embl", "embl"):
          print(record.id)
          print(record.name)
          print(record.description)
        
        # output
        CAA39891.1
        CAA39891
        Lactococcus lactis hypothetical protein

You can also convert files from one format to another:

        from Bio import SeqIO
        count = SeqIO.convert("someEmblFile.gb", "genbank", "someEmblFile.embl", "embl')
        print("Successfully converted %i files" % count)

Source: https://biopython.org/wiki/SeqIO

Alternatively, this markdown file has three parts: [part 1](https://github.com/disulfidebond/Genomic_annotations/blob/master/media/parse_EMBL_file_pt1.md), [part 2](https://github.com/disulfidebond/Genomic_annotations/blob/master/media/parse_EMBL_file_pt2.md), [part 3](https://github.com/disulfidebond/Genomic_annotations/blob/master/media/parse_EMBL_file_pt3.md). It explains how to parse out data from an EMBL file. It is well-commented, and points out where modifications can be made.
