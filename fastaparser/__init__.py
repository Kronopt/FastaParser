#!python
# coding: utf-8

"""
FastaParser

Parses FASTA files with the Reader class (generates FastaSequence objects or plane strings of the parsed sequences)
Writes FASTA files with the Writer class (takes FastaSequence objects or headers + sequences as strings)

ex:
    > import fastaparser
    > with open('fasta_file.fasta') as fasta_file:
    >   reader = fastaparser.Reader(fasta_file)
    >   [seq.id for seq in reader]
    ['HSBGPG', 'HSGLTH1']

ex:
    > import fastaparser
    > with open('fasta_file.fasta', 'w') as fasta_file:
    >   writer = fastaparser.Writer(fasta_file)
    >   seqs = [('HSBGPG example sequence', 'TTCCAGGTGTGCCAATCCAGTCCATG'),
    >   ...     ('HSGLTH1 example sequence 2', 'GTACCTGACCTAACCGTGTGGACCTT')]
    > writer.writefastas(seqs)

Based on these pages:
http://genetics.bwh.harvard.edu/pph/FASTA.html
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
https://en.wikipedia.org/wiki/FASTA_format
https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta
http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
"""

__author__ = "Pedro HC David, https://github.com/Kronopt"
__credits__ = ["Pedro HC David"]
__version__ = "1.1.1"
__license__ = "GPLv3"


from .constants import *
from .fastasequence import FastaSequence
from .lettercode import LetterCode
from .parsedefinitionline import ParseDefinitionLine
from .reader import Reader
from .writer import Writer
