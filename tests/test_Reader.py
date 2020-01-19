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


@pytest.fixture()
def fasta_empty():
    f = open('tests/fasta_empty.fasta', 'r')
    yield f
    f.close()


@pytest.fixture()
def fasta_nucleotide_single():
    f = open('tests/fasta_nucleotide_single.fasta', 'r')
    yield f
    f.close()


@pytest.fixture()
def fasta_nucleotide_multiple():
    f = open('tests/fasta_nucleotide_multiple.fasta', 'r')
    yield f
    f.close()


@pytest.fixture()
def fasta_aminoacid_single():
    f = open('tests/fasta_aminoacid_single.fasta', 'r')
    yield f
    f.close()


@pytest.fixture()
def fasta_aminoacid_multiple():
    f = open('tests/fasta_aminoacid_multiple.fasta', 'r')
    yield f
    f.close()


#######
# Tests
#######


class Test__init__:
    def test_fasta_file_object_good(self, fasta_nucleotide_multiple):
        fasta_reader = Reader(fasta_nucleotide_multiple)
        assert fasta_reader.fasta_file is fasta_nucleotide_multiple
        assert fasta_reader.sequences_type is None
        assert fasta_reader.infer_type is False
        assert fasta_reader.parse_method == 'rich'

    def test_fasta_file_object_closed(self, fasta_nucleotide_multiple):
        fasta_nucleotide_multiple.close()
        with pytest.raises(TypeError):
            Reader(fasta_nucleotide_multiple)

    def test_fasta_file_object_not_a_file(self):
        with pytest.raises(TypeError):
            Reader('')
        with pytest.raises(TypeError):
            Reader([])
        with pytest.raises(TypeError):
            Reader(123)

    def test_sequences_type_nucleotide(self, fasta_nucleotide_multiple):
        fasta_reader = Reader(fasta_nucleotide_multiple, sequences_type='nucleotide')
        assert fasta_reader.fasta_file is fasta_nucleotide_multiple
        assert fasta_reader.sequences_type == 'nucleotide'
        assert fasta_reader.infer_type is False
        assert fasta_reader.parse_method == 'rich'

    def test_sequences_type_aminoacid(self, fasta_aminoacid_multiple):
        fasta_reader = Reader(fasta_aminoacid_multiple, sequences_type='aminoacid')
        assert fasta_reader.fasta_file is fasta_aminoacid_multiple
        assert fasta_reader.sequences_type == 'aminoacid'
        assert fasta_reader.infer_type is False
        assert fasta_reader.parse_method == 'rich'

    # test_sequences_type_non (already tested)

    def test_sequences_type_wrong_type(self, fasta_empty):
        with pytest.raises(TypeError):
            Reader(fasta_empty, sequences_type=[])
        with pytest.raises(TypeError):
            Reader(fasta_empty, sequences_type=123)

    def test_sequences_type_wrong_str(self, fasta_empty):
        with pytest.raises(TypeError):
            Reader(fasta_empty, sequences_type='')
        with pytest.raises(TypeError):
            Reader(fasta_empty, sequences_type=' ')
        with pytest.raises(TypeError):
            Reader(fasta_empty, sequences_type='wrong_type')

    def test_infer_type_true(self, fasta_empty):
        fasta_reader = Reader(fasta_empty, infer_type=True)
        assert fasta_reader.fasta_file is fasta_empty
        assert fasta_reader.sequences_type is None
        assert fasta_reader.infer_type is True
        assert fasta_reader.parse_method == 'rich'

    # test_infer_type_false (already tested)

    def test_infer_type_not_bool(self, fasta_empty):
        with pytest.raises(TypeError):
            Reader(fasta_empty, infer_type='')
        with pytest.raises(TypeError):
            Reader(fasta_empty, infer_type=[])
        with pytest.raises(TypeError):
            Reader(fasta_empty, infer_type=123)

    # test_parse_method_rich (already tested)

    def test_parse_method_quick(self, fasta_empty):
        fasta_reader = Reader(fasta_empty, parse_method='quick')
        assert fasta_reader.fasta_file is fasta_empty
        assert fasta_reader.sequences_type is None
        assert fasta_reader.infer_type is False
        assert fasta_reader.parse_method == 'quick'

    def test_parse_method_wrong_type(self, fasta_empty):
        with pytest.raises(TypeError):
            Reader(fasta_empty, parse_method=[])
        with pytest.raises(TypeError):
            Reader(fasta_empty, parse_method=123)

    def test_parse_method_wrong_str(self, fasta_empty):
        with pytest.raises(TypeError):
            Reader(fasta_empty, parse_method='')
        with pytest.raises(TypeError):
            Reader(fasta_empty, parse_method=' ')
        with pytest.raises(TypeError):
            Reader(fasta_empty, parse_method='wrong_type')


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
