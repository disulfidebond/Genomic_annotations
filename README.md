# Genomic annotations and file formats
A recent project involved manipulating and reformating common Bioinformatics file formats, such as Genbank, EMBL, and fasta. The required tasks were only vaguely connected, so this repository is a cookbook for file manipulation of common Bioinformatics file formats, including parsing and I/O. The general format is to explain the file format, explain possible tasks, and then detail how BioPython can be used to accomplish these goals, and at least one possibility for completing the task from scratch.

## Fasta
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
