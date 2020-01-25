# FastaParser (WORK IN PROGRESS)
A python FASTA parser

#### Install Dependencies
* Python 3.8 (probably 2.7-3.8)

#### How to run
* Download this repo

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
* Documentation
    * Installation
    * Usage
    * Examples
    * API specification
    * Contributing
    * Authors
    * History
* Conda package (?)

#### Maybe
* Identify FASTA ID's
* FASTAID class, to then return in the id property of FastaSequence
* FASTQ parser
* Per fasta sequence, show warning if there are characters not in the FASTA specification
* Way to disable warnings
* Allow setting of id/description
