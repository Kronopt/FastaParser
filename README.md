# PyFastaParser
A simple python FASTA parser

#### Install Dependencies
* Python 3.7 (probably 2.7-3.7)

#### How to run
* Download this repo
* Install dependencies

```Python
>>> import FastaParser
>>> parser = FastaParser.FastaParser("fasta_file.fasta")
>>> for seq in parser:
...     print(seq.id, seq.description)

HSBGPG Human gene for bone gla protein (BGP)
HSGLTH1 Human theta 1-globin gene
```

#### Work in Progress
* Documentation
* Examples
* Tests
* PyPi package
