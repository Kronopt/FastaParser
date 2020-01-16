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
    # test_fasta_file_object_good
    # test_fasta_file_object_closed
    # test_fasta_file_object_not_a_file

    # test_sequences_type_nucleotide
    # test_sequences_type_aminoacid
    # test_sequences_type_none
    # test_sequences_type_wrong_type
    # test_sequences_type_wrong_str

    # test_infer_type_true
    # test_infer_type_false
    # test_infer_type_not_bool

    # test_parse_method_rich
    # test_parse_method_quick
    # test_parse_method_wrong_type
    # test_parse_method_wrong_str
    pass


class Test__iter__:
    pass


class Test__next__:
    pass


class Test__repr__:
    pass


# tested in Test__Init__:
#   class Test_fasta_file_property
#   class Test_sequences_type
#   class Test_infer_type_property
#   class Test_parse_method_property


# TODO
# empty file
# single fasta file
# multiple fasta file

# quick and rich parsing methods
# nucleotide and aminoacid sequences

# empty lines at the end of file
# empty lines at the beginning of file
# no empty lines at the end of file
# empty lines between fastas in multiple fasta file
# empty lines in the middle of a sequence
# headers one after the other
# characters not in fasta specification
