#!python
# coding: utf-8

"""
Tests for fastaparser.Reader class.
"""


import os
import pytest
from fastaparser import Reader
from conftest import fasta_contents


##########
# Fixtures
##########


@pytest.fixture()
def fasta_empty():
    f = open('tests/fasta_empty.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_no_empty_lines():
    f = open('tests/fasta_no_empty_lines.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_no_empty_lines_contents():
    return fasta_contents('tests/fasta_no_empty_lines.fasta')


@pytest.fixture()
def fasta_multiple_empty_lines():
    f = open('tests/fasta_multiple_empty_lines.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_multiple_empty_lines_contents():
    return fasta_contents('tests/fasta_multiple_empty_lines.fasta')


@pytest.fixture()
def fasta_empty_lines_in_sequence():
    f = open('tests/fasta_empty_lines_in_sequence.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_empty_lines_in_sequence_contents():
    return fasta_contents('tests/fasta_empty_lines_in_sequence.fasta')


@pytest.fixture()
def fasta_aminoacid_single():
    f = open('tests/fasta_aminoacid_single.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_aminoacid_single_contents():
    return fasta_contents('tests/fasta_aminoacid_single.fasta')


@pytest.fixture()
def fasta_aminoacid_multiple():
    f = open('tests/fasta_aminoacid_multiple.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_aminoacid_multiple_contents():
    return fasta_contents('tests/fasta_aminoacid_multiple.fasta')


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

    def test_current_iterator(self, fasta_empty):
        fasta_reader = Reader(fasta_empty)
        assert fasta_reader._current_iterator is None


class Test__iter__:
    def test_closed_file(self, fasta_empty):
        fasta_reader = Reader(fasta_empty)
        fasta_empty.close()
        with pytest.raises(TypeError):
            fasta_reader.__iter__()

    def test_current_iterator(self, fasta_empty):
        fasta_reader = Reader(fasta_empty)
        assert fasta_reader._current_iterator is None
        iterator = fasta_reader.__iter__()
        assert fasta_reader._current_iterator == iterator
        try:  # test if is iterable
            iter(fasta_reader._current_iterator)
        except TypeError:
            pytest.fail('Reader._current_iterator is not iterable.')

    def test_empty_fasta_file_rich(self, fasta_empty):
        fasta_reader = Reader(fasta_empty)
        fastas = []
        for fasta in fasta_reader:
            fastas.append(fasta)
        assert len(fastas) == 0

    def test_empty_fasta_file_quick(self, fasta_empty):
        fasta_reader = Reader(fasta_empty, parse_method='quick')
        fastas = []
        for fasta in fasta_reader:
            fastas.append(fasta)
        assert len(fastas) == 0

    def test_single_fasta_file_rich(self, fasta_nucleotide_single, fasta_nucleotide_single_contents,
                                    fasta_aminoacid_single, fasta_aminoacid_single_contents):
        # nucleotide
        fasta_nucleotide_reader = Reader(fasta_nucleotide_single, sequences_type='nucleotide')
        fastas_nucleotide = []
        for fasta in fasta_nucleotide_reader:
            fastas_nucleotide.append(fasta)
        assert len(fastas_nucleotide) == 1
        assert fastas_nucleotide[0].sequence_as_string() == fasta_nucleotide_single_contents[0][2]
        assert fastas_nucleotide[0].id == fasta_nucleotide_single_contents[0][0]
        assert fastas_nucleotide[0].description == fasta_nucleotide_single_contents[0][1]
        assert fastas_nucleotide[0].sequence_type == 'nucleotide'
        assert fastas_nucleotide[0].inferred_type is False
        # aminoacid
        fasta_aminoacid_reader = Reader(fasta_aminoacid_single, sequences_type='aminoacid')
        fastas_aminoacid = []
        for fasta in fasta_aminoacid_reader:
            fastas_aminoacid.append(fasta)
        assert len(fastas_aminoacid) == 1
        assert fastas_aminoacid[0].sequence_as_string() == fasta_aminoacid_single_contents[0][2]
        assert fastas_aminoacid[0].id == fasta_aminoacid_single_contents[0][0]
        assert fastas_aminoacid[0].description == fasta_aminoacid_single_contents[0][1]
        assert fastas_aminoacid[0].sequence_type == 'aminoacid'
        assert fastas_aminoacid[0].inferred_type is False

    def test_single_fasta_file_quick(self, fasta_nucleotide_single, fasta_nucleotide_single_contents,
                                     fasta_aminoacid_single, fasta_aminoacid_single_contents):
        # nucleotide
        fasta_nucleotide_reader = Reader(fasta_nucleotide_single, parse_method='quick')
        fastas_nucleotide = []
        for fasta in fasta_nucleotide_reader:
            fastas_nucleotide.append(fasta)
        assert len(fastas_nucleotide) == 1
        assert fastas_nucleotide[0].sequence == fasta_nucleotide_single_contents[0][2]
        assert fastas_nucleotide[0].header == '>' + ' '.join((fasta_nucleotide_single_contents[0][0],
                                                              fasta_nucleotide_single_contents[0][1]))
        # aminoacid
        fasta_aminoacid_reader = Reader(fasta_aminoacid_single, parse_method='quick')
        fastas_aminoacid = []
        for fasta in fasta_aminoacid_reader:
            fastas_aminoacid.append(fasta)
        assert len(fastas_aminoacid) == 1
        assert fastas_aminoacid[0].sequence == fasta_aminoacid_single_contents[0][2]
        assert fastas_aminoacid[0].header == '>' + ' '.join((fasta_aminoacid_single_contents[0][0],
                                                             fasta_aminoacid_single_contents[0][1]))

    def test_multiple_fasta_file_rich(self, fasta_nucleotide_multiple, fasta_nucleotide_multiple_contents,
                                      fasta_aminoacid_multiple, fasta_aminoacid_multiple_contents):
        # nucleotide
        fasta_nucleotide_reader = Reader(fasta_nucleotide_multiple, sequences_type='nucleotide')
        fastas_nucleotide = []
        for fasta in fasta_nucleotide_reader:
            fastas_nucleotide.append(fasta)
            assert fasta.sequence_as_string() == fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][2]
            assert fasta.id == fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][0]
            assert fasta.description == fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][1]
            assert fasta.sequence_type == 'nucleotide'
            assert fasta.inferred_type is False
        assert len(fastas_nucleotide) == 17
        # aminoacid
        fasta_aminoacid_reader = Reader(fasta_aminoacid_multiple, sequences_type='aminoacid')
        fastas_aminoacid = []
        for fasta in fasta_aminoacid_reader:
            fastas_aminoacid.append(fasta)
            assert fasta.sequence_as_string() == fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][2]
            assert fasta.id == fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][0]
            assert fasta.description == fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][1]
            assert fasta.sequence_type == 'aminoacid'
            assert fasta.inferred_type is False
        assert len(fastas_aminoacid) == 20

    def test_multiple_fasta_file_quick(self, fasta_nucleotide_multiple, fasta_nucleotide_multiple_contents,
                                       fasta_aminoacid_multiple, fasta_aminoacid_multiple_contents):
        # nucleotide
        fasta_nucleotide_reader = Reader(fasta_nucleotide_multiple, parse_method='quick')
        fastas_nucleotide = []
        for fasta in fasta_nucleotide_reader:
            fastas_nucleotide.append(fasta)
            assert fasta.sequence == fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][2]
            assert fasta.header == '>' + ' '.join((fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][0],
                                                   fasta_nucleotide_multiple_contents[len(fastas_nucleotide)-1][1]))
        assert len(fastas_nucleotide) == 17
        # aminoacid
        fasta_aminoacid_reader = Reader(fasta_aminoacid_multiple, parse_method='quick')
        fastas_aminoacid = []
        for fasta in fasta_aminoacid_reader:
            fastas_aminoacid.append(fasta)
            assert fasta.sequence == fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][2]
            assert fasta.header == '>' + ' '.join((fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][0],
                                                   fasta_aminoacid_multiple_contents[len(fastas_aminoacid)-1][1]))
        assert len(fastas_aminoacid) == 20

    def test_empty_lines_between_fastas(self, fasta_multiple_empty_lines, fasta_multiple_empty_lines_contents):
        fasta_reader = Reader(fasta_multiple_empty_lines, parse_method='quick')
        fastas = []
        for fasta in fasta_reader:
            fastas.append(fasta)
            assert fasta.sequence == fasta_multiple_empty_lines_contents[len(fastas)-1][2]
            assert fasta.header == '>' + ' '.join((fasta_multiple_empty_lines_contents[len(fastas)-1][0],
                                                   fasta_multiple_empty_lines_contents[len(fastas)-1][1]))
        assert len(fastas) == 2

    def test_fasta_no_empty_lines(self, fasta_no_empty_lines, fasta_no_empty_lines_contents):
        fasta_reader = Reader(fasta_no_empty_lines, parse_method='quick')
        fastas = []
        for fasta in fasta_reader:
            fastas.append(fasta)
            assert fasta.sequence == fasta_no_empty_lines_contents[len(fastas)-1][2]
            assert fasta.header == '>' + ' '.join((fasta_no_empty_lines_contents[len(fastas)-1][0],
                                                   fasta_no_empty_lines_contents[len(fastas)-1][1]))
        assert len(fastas) == 2

    def test_fasta_empty_lines_in_sequence(self, fasta_empty_lines_in_sequence, fasta_empty_lines_in_sequence_contents):
        fasta_reader = Reader(fasta_empty_lines_in_sequence, parse_method='quick')
        fastas = []
        for fasta in fasta_reader:
            fastas.append(fasta)
            assert fasta.sequence == fasta_empty_lines_in_sequence_contents[len(fastas)-1][2]
            assert fasta.header == '>' + ' '.join((fasta_empty_lines_in_sequence_contents[len(fastas)-1][0],
                                                   fasta_empty_lines_in_sequence_contents[len(fastas)-1][1]))
        assert len(fastas) == 2


class Test__next__:
    def test_existing_current_iterator(self, fasta_nucleotide_multiple, fasta_nucleotide_multiple_contents):
        fasta_reader = Reader(fasta_nucleotide_multiple, sequences_type='nucleotide')
        iter(fasta_reader)
        fasta = next(fasta_reader)
        assert fasta.sequence_as_string() == fasta_nucleotide_multiple_contents[0][2]
        assert fasta.id == fasta_nucleotide_multiple_contents[0][0]
        assert fasta.description == fasta_nucleotide_multiple_contents[0][1]
        assert fasta.sequence_type == 'nucleotide'
        assert fasta.inferred_type is False
        fasta = next(fasta_reader)
        assert fasta.sequence_as_string() == fasta_nucleotide_multiple_contents[1][2]
        assert fasta.id == fasta_nucleotide_multiple_contents[1][0]
        assert fasta.description == fasta_nucleotide_multiple_contents[1][1]
        assert fasta.sequence_type == 'nucleotide'
        assert fasta.inferred_type is False

    def test_no_current_iterator(self, fasta_nucleotide_multiple, fasta_nucleotide_multiple_contents):
        fasta_reader = Reader(fasta_nucleotide_multiple, sequences_type='nucleotide')
        fasta = next(fasta_reader)
        assert fasta.sequence_as_string() == fasta_nucleotide_multiple_contents[0][2]
        assert fasta.id == fasta_nucleotide_multiple_contents[0][0]
        assert fasta.description == fasta_nucleotide_multiple_contents[0][1]
        assert fasta.sequence_type == 'nucleotide'
        assert fasta.inferred_type is False
        fasta = next(fasta_reader)
        assert fasta.sequence_as_string() == fasta_nucleotide_multiple_contents[1][2]
        assert fasta.id == fasta_nucleotide_multiple_contents[1][0]
        assert fasta.description == fasta_nucleotide_multiple_contents[1][1]
        assert fasta.sequence_type == 'nucleotide'
        assert fasta.inferred_type is False


class Test__repr__:
    def test__repr__(self, fasta_empty):
        fasta_reader = Reader(fasta_empty)
        assert repr(fasta_reader) == 'fastaparser.Reader(%s)' % os.path.abspath(fasta_empty.name)


# tested in Test__Init__:
#   class Test_fasta_file_property
#   class Test_sequences_type
#   class Test_infer_type_property
#   class Test_parse_method_property
