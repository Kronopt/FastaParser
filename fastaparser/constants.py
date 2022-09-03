#!python
# coding: utf-8

"""
Constants used throughout the whole package.
These include letter codes and their respective names for all "good" and degenerate nucleotides/aminoacids.
"""

# NUCLEOTIDE DICTIONARIES
NUCLEOTIDE_LETTER_CODES_GOOD = {
    "A": "adenosine",
    "C": "cytidine",
    "G": "guanine",
    "T": "thymidine",
    "N": "any (A/G/C/T)",
    "U": "uridine",
}
NUCLEOTIDE_LETTER_CODES_DEGENERATE = {
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
NUCLEOTIDE_LETTER_CODES_COMPLEMENT = {
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

# AMINOACID DICTIONARIES
AMINOACID_LETTER_CODES_GOOD = {
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
AMINOACID_LETTER_CODES_DEGENERATE = {"-": "gap of indeterminate length"}

# BOTH NUCLEOTIDE AND AMINOACID
LETTER_CODES = {
    "nucleotide": (NUCLEOTIDE_LETTER_CODES_GOOD, NUCLEOTIDE_LETTER_CODES_DEGENERATE),
    "aminoacid": (AMINOACID_LETTER_CODES_GOOD, AMINOACID_LETTER_CODES_DEGENERATE),
}

# set operations
LETTER_CODES_ALL = set(
    list(NUCLEOTIDE_LETTER_CODES_GOOD)
    + list(NUCLEOTIDE_LETTER_CODES_DEGENERATE)
    + list(AMINOACID_LETTER_CODES_GOOD)
    + list(AMINOACID_LETTER_CODES_DEGENERATE)
)
NUCLEOTIDE_LETTER_CODES_ALL = set(
    list(NUCLEOTIDE_LETTER_CODES_GOOD) + list(NUCLEOTIDE_LETTER_CODES_DEGENERATE)
)
AMINOACID_LETTER_CODES_ALL = set(
    list(AMINOACID_LETTER_CODES_GOOD) + list(AMINOACID_LETTER_CODES_DEGENERATE)
)
AMINOACIDS_NOT_IN_NUCLEOTIDES = AMINOACID_LETTER_CODES_ALL - NUCLEOTIDE_LETTER_CODES_ALL
# nucleotides_not_in_aminoacids would be empty
