#!python
# coding: utf-8

"""
Writer - FASTA writer.
"""

import os
from .fastasequence import FastaSequence
from .parsedefinitionline import ParseDefinitionLine


class Writer(ParseDefinitionLine):
    """
    Writer for the given FASTA file.
    Writes FastaSequence objects or tuples of (header, sequence) to the given file.

    Attributes
    ----------
    fasta_file : file object
        The FASTA file passed as parameter.

    Methods
    -------
    writefasta(FastaSequence or (header, sequence))
        Writes a single FASTA sequence.
        The FASTA sequence can either be a FastaSequence object or a tuple of (header, sequence) as strings (in this
        order).
        header can be an empty string.
        sequence must have content.
    writefastas(list of: FastaSequence or (header, sequence))
        Writes multiple FASTA sequences.
        FASTA sequences in the list should be defined as in writefasta().

    Raises
    ------
    TypeError
        When calling __init__, if fasta_file is of the wrong type.
        When calling __init__, if fasta_file is not a file object, is closed or is not writable.
        When calling writefasta(), if fasta_sequence is of the wrong type.
        When calling writefastas(), if fasta_sequences is not iterable.
    """

    def __init__(self, fasta_file):
        """
        Initializes file object (checks if fasta_file is a file object opened for writing).

        Parameters
        ----------
        fasta_file : file object
            An opened file handle ready for writing.

        Raises
        ------
        TypeError
            If fasta_file is of the wrong type.
            If fasta_file is not a file object, is closed or is not writable.
        """
        # assume it's a file object
        if (
            hasattr(fasta_file, "writelines")
            and hasattr(fasta_file, "closed")
            and hasattr(fasta_file, "writable")
        ):
            if not fasta_file.closed and fasta_file.writable():
                self._fasta_file = fasta_file
            else:
                raise TypeError("fasta_file must be opened for writing")
        else:
            raise TypeError("fasta_file must be a file object")

    @property
    def fasta_file(self):
        """return fasta_file."""
        return self._fasta_file

    def writefasta(self, fasta_sequence):
        """
        Writes a single FASTA sequence to the provided file.
        Open file with mode 'a' to append sequences to an existing FASTA file.

        Parameters
        ----------
        fasta_sequence : FastaSequence or (header : str, sequence : str)
            A FASTA sequence is built from the data contained in the provided FastaSequence object or the tuple of
            header + sequence.
            header may contain or not the starting '>'. header can be an empty string.
            sequence must be a non empty string.

        Raises
        ------
        TypeError
            If fasta_sequence is of the wrong type.
        """
        # either use the FastaSequence object directly
        if isinstance(fasta_sequence, FastaSequence):
            pass

        # or create one with the provided header and sequence
        elif (
            isinstance(fasta_sequence, (tuple, list))
            and len(fasta_sequence) == 2
            and isinstance(fasta_sequence[0], str)
            and isinstance(fasta_sequence[1], str)
        ):
            id_, description = self._parse_definition_line(fasta_sequence[0])
            sequence = "".join(
                fasta_sequence[1].split("\n")
            )  # remove '\n's from sequence
            fasta_sequence = FastaSequence(sequence, id_, description)

        else:
            raise TypeError(
                "fasta_sequence must be a FastaSequence object or a tuple (header : str, sequence : str)"
            )

        # write fasta to file
        self._fasta_file.write(fasta_sequence.formatted_fasta() + "\n\n")

    def writefastas(self, fasta_sequences):
        """
        Writes multiple FASTA sequences to the provided file.
        Simply calls writefasta() for each object in fasta_sequences.
        Open the file with mode 'a' if you want to append multiple sequences to an existing FASTA file.

        Parameters
        ----------
        fasta_sequences : iterable of FastaSequence or iterable of (header : str, sequence : str)
            FASTA sequences are built from the data contained in the provided FastaSequence objects or the tuples of
            header + sequence.
            headers may contain or not the starting '>'. headers can be empty strings.
            sequences must be non empty strings.

        Raises
        ------
        TypeError
            If fasta_sequences is not iterable.
        """
        try:
            iter(fasta_sequences)
        except TypeError as exc:
            raise TypeError(
                "fasta_sequences must be an iterable of FastaSequence "
                "objects or an iterable of tuples (header : str, sequence : str)"
            ) from exc

        for fasta in fasta_sequences:
            self.writefasta(fasta)

    def __repr__(self):
        return "fastaparser.Writer(%s)" % os.path.abspath(self._fasta_file.name)
