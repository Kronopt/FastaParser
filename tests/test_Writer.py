#!python
# coding: utf-8

"""
Tests for FastaParser.Writer class
"""


import hashlib
import os
import pytest
from FastaParser import Reader, Writer


##########
# Fixtures
##########


def compare_2_files(file_read, file_written):
    hash_read = hashlib.md5()
    file_read.seek(0)  # reset cursor of read file
    for line in file_read:
        hash_read.update(line.encode())

    hash_written = hashlib.md5()
    file_written.close()  # flush buffer to temporary file
    with open(file_written.name) as file_written_read:  # re-open written file as read
        for line in file_written_read:
            hash_written.update(line.encode())

    assert hash_read.hexdigest() == hash_written.hexdigest()


@pytest.fixture()
def fasta_temporary_file():
    f = open('tests/FASTA_TEMPORARY_FILE_FOR_WRITING.fasta', 'w')
    yield f
    f.close()
    os.remove(f.name)


#######
# Tests
#######


class Test__init__:
    def test_fasta_file_object_good(self, fasta_temporary_file):
        fasta_writer = Writer(fasta_temporary_file)
        assert fasta_writer.fasta_file is fasta_temporary_file

    def test_fasta_file_object_closed(self, fasta_temporary_file):
        fasta_temporary_file.close()
        with pytest.raises(TypeError):
            Writer(fasta_temporary_file)

    def test_fasta_file_object_not_a_file(self):
        with pytest.raises(TypeError):
            Writer('')
        with pytest.raises(TypeError):
            Writer([])
        with pytest.raises(TypeError):
            Writer(123)


class Test_writefasta:
    def test_fasta_sequence_fastasequence_object(self, fasta_nucleotide_single, fasta_temporary_file):
        fasta_reader = Reader(fasta_nucleotide_single)
        fasta_writer = Writer(fasta_temporary_file)
        fasta_writer.writefasta(next(fasta_reader))  # only a single FASTA sequence in fasta_nucleotide_single
        # at this point the 2 files should be equal
        compare_2_files(fasta_nucleotide_single, fasta_temporary_file)

    def test_fasta_sequence_tuple(self, fasta_nucleotide_single, fasta_temporary_file):
        fasta_reader = Reader(fasta_nucleotide_single)
        fasta_writer = Writer(fasta_temporary_file)
        fasta = next(fasta_reader)
        fasta_writer.writefasta((fasta.formatted_definition_line(), fasta.formatted_sequence()))
        # at this point the 2 files should be equal
        compare_2_files(fasta_nucleotide_single, fasta_temporary_file)

    def test_fasta_sequence_wrong_type(self, fasta_temporary_file):
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefasta('')
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefasta(123)
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefasta([])
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefasta((1, 2))


class Test_writefastas:
    def test_fasta_sequence_fastasequence_objects(self, fasta_nucleotide_multiple, fasta_temporary_file):
        fasta_reader = Reader(fasta_nucleotide_multiple)
        fasta_writer = Writer(fasta_temporary_file)
        fasta_writer.writefastas(fasta_reader)
        # at this point the 2 files should be equal
        compare_2_files(fasta_nucleotide_multiple, fasta_temporary_file)

    def test_fasta_sequence_tuples(self, fasta_nucleotide_multiple, fasta_temporary_file):
        fasta_reader = Reader(fasta_nucleotide_multiple)
        fasta_writer = Writer(fasta_temporary_file)
        fastas = [(fasta.formatted_definition_line(), fasta.formatted_sequence()) for fasta in fasta_reader]
        fasta_writer.writefastas(fastas)
        # at this point the 2 files should be equal
        compare_2_files(fasta_nucleotide_multiple, fasta_temporary_file)

    def test_fasta_sequence_wrong_type(self, fasta_temporary_file):
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefastas(123)
            with pytest.raises(TypeError):
                fasta_writer = Writer(fasta_temporary_file)
                fasta_writer.writefastas(None)
        with pytest.raises(TypeError):
            fasta_writer = Writer(fasta_temporary_file)
            fasta_writer.writefastas([1, 2])


class Test__repr__:
    def test__repr__(self, fasta_temporary_file):
        fasta_writer = Writer(fasta_temporary_file)
        assert repr(fasta_writer) == 'FastaParser.FastaParser.Writer(%s)' % os.path.abspath(fasta_temporary_file.name)


# tested in Test__Init__:
#   class Test_fasta_file_property
