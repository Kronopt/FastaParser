# PyFastaParser
A simple python FASTA parser (WORK IN PROGRESS)

#### Install Dependencies
* Python 3.7 (probably 2.7-3.7)

#### How to run
* Download this repo

```Python
>>> import FastaParser
>>> with open("fasta_file.fasta") as fasta_file:
...     parser = FastaParser.Reader(fasta_file)
...     for seq in parser:
...         print(seq.id, seq.description)

HSBGPG Human gene for bone gla protein (BGP)
HSGLTH1 Human theta 1-globin gene
```

#### Work in Progress
* Asyncio for Reader/Writer
* Context Manager (?)
* Documentation
* Examples
* Read the whole fasta "specification" to provide more features
* Tests
* PyPi package
