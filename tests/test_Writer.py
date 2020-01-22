#!python
# coding: utf-8

"""
Tests for FastaParser.Writer class
"""


import os
import pytest
from FastaParser import Writer


##########
# Fixtures
##########


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


# TODO
# class Test_writefasta:
    # write FastaSequence
    # write (header, sequence)
# class Test_writefastas:
# class Test__repr__:




# tested in Test__Init__:
#   class Test_fasta_file_property
