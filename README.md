# FastaParser
A python FASTA parser

```Python
>>> import fastaparser
>>> with open("fasta_file.fasta") as fasta_file:
...     parser = fastaparser.Reader(fasta_file)
...     for seq in parser:
...         print(seq.id, seq.description)

HSBGPG Human gene for bone gla protein (BGP)
HSGLTH1 Human theta 1-globin gene
```

#### Work in Progress
* PyPi package
* Documentation (readthedocs)
    * Home
    * Installation
    * Usage
    * Examples
    * API Specification
    * Contributing
    * Authors
    * History
* Conda package (?)
* README

#### Maybe
* Identify FASTA ID's
* FASTAID class, to then return in the id property of FastaSequence
* FASTQ parser
* Per fasta sequence, show warning if there are characters not in the FASTA specification
* Way to disable warnings
* Allow setting of id/description
