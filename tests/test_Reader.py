#!python
# coding: utf-8

"""
Tests for FastaParser.Reader class
"""


import pytest
from FastaParser import Reader


##########
# Fixtures
##########


#######
# Tests
#######


class Test__init__:
    pass


# TODO
# empty file
# single fasta file
# multiple fasta file
# empty lines at the end of file
# empty lines at the beginning of file
# no empty lines at the end of file
# empty lines between fastas in multiple fasta file
# empty lines in the middle of a sequence
# headers one after the other
# characters not in fasta specification
