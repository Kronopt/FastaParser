# Examples
Assume `fasta_file.fasta` exists and contains 2 FASTA sequences.

## Read FASTA files
Generate python objects from FASTA files:

```Python
>>> import fastaparser
>>> with open("fasta_file.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            print('ID:', seq.id)
            print('Description:', seq.description)
            print('Sequence:', seq.sequence_as_string())
            print()
```
output:
```
ID: sp|P04439|HLAA_HUMAN
Description: HLA class I histocompatibility antigen, A alpha chain OS=Homo sapi...
Sequence: MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRM...

ID: sp|P15822|ZEP1_HUMAN
Description: Zinc finger protein 40 OS=Homo sapiens OX=9606 GN=HIVEP1 PE=1 SV=3...
Sequence: MPRTKQIHPRNLRDKIEEAQKELNGAEVSKKEILQAGVKGTSESLKGVKRKKIVAENHLKKIPKSPLRN...
```

or just parse FASTA headers and sequences, which is much faster but less feature rich:
```Python
>>> import fastaparser
>>> with open("fasta_file.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file, parse_method='quick')
        for seq in parser:
            # seq is a namedtuple('Fasta', ['header', 'sequence'])
            print('Header:', seq.header)
            print('Sequence:', seq.sequence)
            print()
```
output:
```
Header: >sp|P04439|HLAA_HUMAN HLA class I histocompatibility antigen, A alpha c...
Sequence: MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRM...

Header: >sp|P15822|ZEP1_HUMAN Zinc finger protein 40 OS=Homo sapiens OX=9606 GN...
Sequence: MPRTKQIHPRNLRDKIEEAQKELNGAEVSKKEILQAGVKGTSESLKGVKRKKIVAENHLKKIPKSPLRN...
```

## Write FASTA files
Create FASTA files from FastaSequence objects:
```Python
>>> import fastaparser
>>> with open("fasta_file.fasta", 'w') as fasta_file:
        writer = fastaparser.Writer(fasta_file)
        fasta_sequence = fastaparser.FastaSequence(
            sequence='ACTGCTGCTAGCTAGC',
            id='id123',
            description='test sequence'
        )
        writer.writefasta(fasta_sequence)
```
or single header and sequence strings:
```Python
>>> import fastaparser
>>> with open("fasta_file.fasta", 'w') as fasta_file:
        writer = fastaparser.Writer(fasta_file)
        writer.writefasta(('id123 test sequence', 'ACTGCTGCTAGCTAGC'))
```
