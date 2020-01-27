# FastaParser

[![pypi](https://img.shields.io/pypi/v/fastaparser "pypi package")](https://pypi.org/project/fastaparser)
[![python versions](https://img.shields.io/pypi/pyversions/fastaparser "supported python versions")](https://pypi.org/project/fastaparser)
[![downloads](https://img.shields.io/pypi/dm/fastaparser "pypi downloads")](https://pypi.org/project/fastaparser)
[![build status](https://github.com/Kronopt/FastaParser/workflows/CI/badge.svg "build status")](https://github.com/Kronopt/FastaParser/actions?query=workflow%3ACI)
[![coverage](https://codecov.io/gh/Kronopt/FastaParser/branch/master/graph/badge.svg "code coverage")](https://codecov.io/gh/Kronopt/FastaParser)
[![license](https://img.shields.io/pypi/l/fastaparser "license")](https://github.com/Kronopt/fastaparser/blob/master/LICENSE)

A python FASTA parser

## Installation
```sh
$ pip install fastaparser
```

## Usage
Generate python objets from FASTA files:
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

Header: >sp|P04439|HLAA_HUMAN HLA class I histocompatibility antigen, A alpha c...
Sequence: MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRM...

Header: >sp|P15822|ZEP1_HUMAN Zinc finger protein 40 OS=Homo sapiens OX=9606 GN...
Sequence: MPRTKQIHPRNLRDKIEEAQKELNGAEVSKKEILQAGVKGTSESLKGVKRKKIVAENHLKKIPKSPLRN...
```

## Documentation
Documentation for FastaParser is available here: [fastaparser.rtfd.io](https://fastaparser.readthedocs.io/en/latest/)

## To do
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

#### Maybe
* Identify FASTA ID's
* FASTAID class, to then return in the id property of FastaSequence
* FASTQ parser
* Per fasta sequence, show warning if there are characters not in the FASTA specification
* Way to disable warnings
* Allow setting of id/description
