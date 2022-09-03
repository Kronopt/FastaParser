#!python
# coding: utf-8

"""
Tests for FastaParser constants.
"""


import pytest
import fastaparser


##########
# Fixtures
##########


@pytest.fixture(scope="module")
def nucleotide_letter_codes_good():
    return {
        "A": "adenosine",
        "C": "cytidine",
        "G": "guanine",
        "T": "thymidine",
        "N": "any (A/G/C/T)",
        "U": "uridine",
    }


@pytest.fixture(scope="module")
def nucleotide_letter_codes_degenerate():
    return {
        "K": "keto (G/T)",
        "S": "strong (G/C)",
        "Y": "pyrimidine (T/C)",
        "M": "amino (A/C)",
        "W": "weak (A/T)",
        "R": "purine (G/A)",
        "B": "G/T/C",
        "D": "G/A/T",
        "H": "A/C/T",
        "V": "G/C/A",
        "-": "gap of indeterminate length",
    }


@pytest.fixture(scope="module")
def nucleotide_letter_codes_complement():
    return {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "U": "A",
        "K": "M",
        "S": "S",
        "Y": "R",
        "M": "K",
        "W": "W",
        "R": "Y",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "-": "-",
    }


@pytest.fixture(scope="module")
def nucleotide_letter_codes_all(
    nucleotide_letter_codes_good, nucleotide_letter_codes_degenerate
):
    return set(
        list(nucleotide_letter_codes_good) + list(nucleotide_letter_codes_degenerate)
    )


@pytest.fixture(scope="module")
def aminoacid_letter_codes_good():
    return {
        "A": "alanine",
        "B": "aspartate/asparagine",
        "C": "cystine",
        "D": "aspartate",
        "E": "glutamate",
        "F": "phenylalanine",
        "G": "glycine",
        "H": "histidine",
        "I": "isoleucine",
        "K": "lysine",
        "L": "leucine",
        "M": "methionine",
        "N": "asparagine",
        "P": "proline",
        "Q": "glutamine",
        "R": "arginine",
        "S": "serine",
        "T": "threonine",
        "U": "selenocysteine",
        "V": "valine",
        "W": "tryptophan",
        "Y": "tyrosine",
        "Z": "glutamate/glutamine",
        "X": "any",
        "*": "translation stop",
    }


@pytest.fixture(scope="module")
def aminoacid_letter_codes_degenerate():
    return {"-": "gap of indeterminate length"}


@pytest.fixture(scope="module")
def aminoacid_letter_codes_all(
    aminoacid_letter_codes_good, aminoacid_letter_codes_degenerate
):
    return set(
        list(aminoacid_letter_codes_good) + list(aminoacid_letter_codes_degenerate)
    )


#######
# Tests
#######


def test_NUCLEOTIDE_LETTER_CODES_GOOD(nucleotide_letter_codes_good):
    assert fastaparser.NUCLEOTIDE_LETTER_CODES_GOOD == nucleotide_letter_codes_good


def test_NUCLEOTIDE_LETTER_CODES_DEGENERATE(nucleotide_letter_codes_degenerate):
    assert (
        fastaparser.NUCLEOTIDE_LETTER_CODES_DEGENERATE
        == nucleotide_letter_codes_degenerate
    )


def test_NUCLEOTIDE_LETTER_CODES_COMPLEMENT(nucleotide_letter_codes_complement):
    assert (
        fastaparser.NUCLEOTIDE_LETTER_CODES_COMPLEMENT
        == nucleotide_letter_codes_complement
    )


def test_AMINOACID_LETTER_CODES_GOOD(aminoacid_letter_codes_good):
    assert fastaparser.AMINOACID_LETTER_CODES_GOOD == aminoacid_letter_codes_good


def test_AMINOACID_LETTER_CODES_DEGENERATE(aminoacid_letter_codes_degenerate):
    assert (
        fastaparser.AMINOACID_LETTER_CODES_DEGENERATE
        == aminoacid_letter_codes_degenerate
    )


def test_LETTER_CODES(
    nucleotide_letter_codes_good,
    nucleotide_letter_codes_degenerate,
    aminoacid_letter_codes_good,
    aminoacid_letter_codes_degenerate,
):
    letter_codes = {
        "nucleotide": (
            nucleotide_letter_codes_good,
            nucleotide_letter_codes_degenerate,
        ),
        "aminoacid": (aminoacid_letter_codes_good, aminoacid_letter_codes_degenerate),
    }
    assert fastaparser.LETTER_CODES == letter_codes


def test_LETTER_CODES_ALL(
    nucleotide_letter_codes_good,
    nucleotide_letter_codes_degenerate,
    aminoacid_letter_codes_good,
    aminoacid_letter_codes_degenerate,
):
    letter_codes_all = set(
        list(nucleotide_letter_codes_good)
        + list(nucleotide_letter_codes_degenerate)
        + list(aminoacid_letter_codes_good)
        + list(aminoacid_letter_codes_degenerate)
    )
    assert fastaparser.LETTER_CODES_ALL == letter_codes_all


def test_NUCLEOTIDE_LETTER_CODES_ALL(nucleotide_letter_codes_all):
    assert fastaparser.NUCLEOTIDE_LETTER_CODES_ALL == nucleotide_letter_codes_all


def test_AMINOACID_LETTER_CODES_ALL(aminoacid_letter_codes_all):
    assert fastaparser.AMINOACID_LETTER_CODES_ALL == aminoacid_letter_codes_all


def test_AMINOACIDS_NOT_IN_NUCLEOTIDES(
    nucleotide_letter_codes_all, aminoacid_letter_codes_all
):
    aminoacids_not_in_nucleotides = (
        aminoacid_letter_codes_all - nucleotide_letter_codes_all
    )
    assert fastaparser.AMINOACIDS_NOT_IN_NUCLEOTIDES == aminoacids_not_in_nucleotides
