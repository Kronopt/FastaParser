#!python
# coding: utf-8

"""
Pytest configuration.
"""


import pytest


def fasta_contents(fasta_path):
    with open(fasta_path) as f:
        fastas = []
        _id = ''
        fasta = ''
        start_of_fasta = True
        for line in f:
            if line.startswith('>') and start_of_fasta:
                start_of_fasta = False
                _id, description = line.rstrip()[1:].split(maxsplit=1)
            elif line.startswith('>'):
                fastas.append((_id, description, fasta))
                _id, description = line.rstrip()[1:].split(maxsplit=1)
                fasta = ''
            else:
                fasta += line.rstrip()
        # last fasta
        fastas.append((_id, description, fasta))
    return fastas  # (id, description, fasta)


@pytest.fixture(scope='session')
def unknown_characters():
    return 'O', '»', '%', 'º', '?', 'غ', '\n'


@pytest.fixture()
def fasta_nucleotide_single():
    f = open('tests/fasta_nucleotide_single.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_nucleotide_single_contents():
    return fasta_contents('tests/fasta_nucleotide_single.fasta')


@pytest.fixture()
def fasta_nucleotide_multiple():
    f = open('tests/fasta_nucleotide_multiple.fasta')
    yield f
    f.close()


@pytest.fixture()
def fasta_nucleotide_multiple_contents():
    return fasta_contents('tests/fasta_nucleotide_multiple.fasta')
