#!python
# coding: utf-8

"""
Reader - FASTA parser/reader.
"""

import os
from collections import namedtuple
from .constants import LETTER_CODES
from .fastasequence import FastaSequence
from .parsedefinitionline import ParseDefinitionLine


class Reader(ParseDefinitionLine):
    """
    Parser/Reader for the given FASTA file.
    Iterates over the FASTA file using one of two parsing mechanisms:
        'rich':
            Returns FastaSequence objects (default).
            Slower, but feature rich.
        'quick':
            Generates objects containing just the FASTA header and sequence attributes
            for each sequence in the FASTA file.
            Parses FASTA files faster but lacks some features.

    Attributes
    ----------
    fasta_file : file object
        The FASTA file passed as parameter.
    sequences_type : 'nucleotide', 'aminoacid' or None
        Indicates the type of sequences to expect ('nucleotide' or 'aminoacid'). Can be None if not known.
    infer_type: bool
        True if Reader was set to infer the sequence type, False otherwise.
    parse_method: 'rich' or 'quick'
        Parse method used ('rich' or 'quick').

    Raises
    ------
    TypeError
        When calling __init__, if fasta_file, sequences_type, infer_type or parse_method are of the wrong type.
        When calling __init__, if fasta_file is not a file object, is closed or is not readable.
        When calling __iter__, if fasta_file is closed.
    """

    _PARSE_METHODS = ("rich", "quick")

    def __init__(
        self, fasta_file, sequences_type=None, infer_type=False, parse_method="rich"
    ):
        """
        Initializes file object (checks if fasta_file is an opened file object).

        Parameters
        ----------
        fasta_file : file object
            An opened file handle for reading.
        sequences_type : 'nucleotide', 'aminoacid' or None, optional
            Indicates the type of sequences to expect ('nucleotide' or 'aminoacid'). None if unknown.
        infer_type : bool, optional
            Indicates if Reader should try to infer aminoacid sequence type for each sequence.
            Can only identify aminoacid sequences.
        parse_method: 'rich' or 'quick', optional
            Parse method to use ('rich' or 'quick'). Defaults to 'rich'.
            'quick' parsing method just parses the header and the sequence into individual properties,
            so it's much faster and less memory intensive. If selected, sequences_type and
            infer_type parameters are ignored.
            'rich' implements more functionality (FastaSequence and LetterCode), but is slower.

        Raises
        ------
        TypeError
            If fasta_file, sequences_type, infer_type or parse_method are of the wrong type.
            If fasta_file is not a file object, is closed or is not readable.
        """
        # for 'quick' parse method
        self._fasta_sequence = namedtuple("Fasta", ["header", "sequence"])

        # assume it's a file object
        if (
            hasattr(fasta_file, "readline")
            and hasattr(fasta_file, "closed")
            and hasattr(fasta_file, "readable")
        ):
            if not fasta_file.closed and fasta_file.readable():
                self._fasta_file = fasta_file
            else:
                raise TypeError("fasta_file must be opened for reading")
        else:
            raise TypeError("fasta_file must be a file object")

        if (
            isinstance(sequences_type, str) and sequences_type in LETTER_CODES
        ) or sequences_type is None:
            self._sequences_type = sequences_type
        else:
            raise TypeError("sequence_type must be one of: %s or None" % LETTER_CODES)

        if isinstance(infer_type, bool):
            self._infer_type = infer_type
        else:
            raise TypeError("infer_type must be bool")

        if isinstance(parse_method, str) and parse_method in self._PARSE_METHODS:
            self._parse_method = parse_method
        else:
            raise TypeError(
                "parse_method must be one of: %s" % ", ".join(self._PARSE_METHODS)
            )

        self._current_iterator = None

    @property
    def fasta_file(self):
        """return fasta_file."""
        return self._fasta_file

    @property
    def sequences_type(self):
        """return sequences_type."""
        return self._sequences_type

    @property
    def infer_type(self):
        """return infer_type."""
        return self._infer_type

    @property
    def parse_method(self):
        """return parse_method."""
        return self._parse_method

    def _generate_fasta_sequence_object(self, sequence, definition_line):
        """
        Generates either a FastaSequence or a namedtuple('Fasta', ['header', 'sequence']) object,
        based on the value of self._parse_method

        Parameters
        ----------
        sequence : str
            Sequence as string.
        definition_line : str
            Definition line (id + description) including '>' at the beginning.

        Returns
        -------
        FastaSequence or namedtuple('Fasta', ['header', 'sequence'])
        """
        if self._parse_method == "rich":
            id_, description = self._parse_definition_line(definition_line)
            fasta_sequence = FastaSequence(
                sequence, id_, description, self._sequences_type, self._infer_type
            )
        else:  # 'quick'
            fasta_sequence = self._fasta_sequence(definition_line, sequence)
        return fasta_sequence

    def _iter_fasta_file(self, fasta_file):
        """
        Iterator of FASTA files (called by __iter__).

        Parameters
        ----------
        fasta_file : file object
            An opened file handle.
        """
        fasta_file.seek(0)  # restart cursor position (just in case)

        definition_line = ""
        sequence = ""

        parsing_fasta_sequence = False
        for line in fasta_file:
            line = line.strip()

            # searching for '>' character at the start of a line
            # following lines contain the sequence, so parsing_fasta_sequence becomes True
            if not parsing_fasta_sequence:
                if line.startswith(">"):
                    definition_line = line
                    parsing_fasta_sequence = True

            # in the middle of parsing a FASTA sequence
            else:
                # builds sequence string until a '>' is found.
                # also ignores empty lines, which is not specified in the FASTA specification,
                # but is forgiven to badly constructed FASTA files
                if len(line) > 0 and line[0] != ">":
                    sequence += line
                elif len(line) > 0 and line[0] == ">":
                    fasta_sequence = self._generate_fasta_sequence_object(
                        sequence, definition_line
                    )
                    yield fasta_sequence

                    # restart variables
                    definition_line = line
                    sequence = ""
                    parsing_fasta_sequence = True

        # end of file, therefore yield last FASTA sequence
        if (
            len(sequence) > 0
        ):  # a FASTA sequence was actually parsed and were not just blank lines
            fasta_sequence = self._generate_fasta_sequence_object(
                sequence, definition_line
            )
            yield fasta_sequence

    def __iter__(self):
        """
        Iterates over the FASTA file.
        Returns a new iterator of the file (from the beginning) every time __iter__ is called.
        """
        if (
            not self._fasta_file.closed and self._fasta_file.readable()
        ):  # check if file is closed
            self._current_iterator = self._iter_fasta_file(self._fasta_file)
            return self._current_iterator
        raise TypeError("fasta_file must be opened for reading")

    def __next__(self):
        """
        Returns the next FASTA sequence from the current iterator (most recent iterator).
        If no iterator still exists, calls __iter__ to create it.
        """
        if self._current_iterator is None:
            self.__iter__()
        return next(self._current_iterator)

    def __repr__(self):
        return "fastaparser.Reader(%s)" % os.path.abspath(self._fasta_file.name)
