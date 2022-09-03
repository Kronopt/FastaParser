#!python
# coding: utf-8

"""
LetterCode - Represents a single nucleotide or aminoacid letter code.
"""

import warnings
from .constants import (
    LETTER_CODES,
    LETTER_CODES_ALL,
    NUCLEOTIDE_LETTER_CODES_COMPLEMENT,
)


warnings.simplefilter(
    "always"
)  # show warnings everytime instead of only the first time they happen


class LetterCode:
    """
    Represents a single letter code.

    Attributes
    ----------
    letter_code : str
        Upper case letter code.
    letter_type : str or None
        'nucleotide' or 'aminoacid'. None if there is no information about sequence type.
    description : str
        Description or nucleotide/aminoacid name of letter code (can be an empty string).
    degenerate : bool or None
        Indicates if a letter code is degenerate or not (can be None if letter code is
        not defined in the FASTA specification or letter_type is unknown).
    supported : bool
        Indicates if letter code is supported or not
        (ie, if letter_type is provided and letter code is defined in the FASTA specification).
    in_fasta_spec : bool
        Indicates if Letter code is defined in the FASTA specification.

    Methods
    -------
    from_lettercode(lettercode)
        Alternate __init__ method. Initializes instance with a LetterCode object as only parameter.
    complement()
        Returns the complementary LetterCode (ideally, of a nucleotide).

    Raises
    ------
    TypeError
        When calling __init__, if letter_code or letter_type are of the wrong type.
        When calling from_lettercode(), if lettercode is of the wrong type.
        When setting letter_type, if letter_type_value is of the wrong type.
        When calling complement(), if letter_type is 'aminoacid'.
    """

    def __init__(self, letter_code, letter_type=None):
        """
        Initializes given letter code.

        Parameters
        ----------
        letter_code : str
            Letter code.
        letter_type : 'nucleotide', 'aminoacid' or None, optional
            'nucleotide' or 'aminoacid' type letter code, None if there is no information.

        Raises
        ------
        TypeError
            If letter_code or letter_type are of the wrong type.
        """
        if isinstance(letter_code, str) and len(letter_code) == 1:
            self._letter_code = letter_code.upper()
            self._in_fasta_spec = self._letter_code in LETTER_CODES_ALL
        else:
            raise TypeError("letter_code must be a single character str")

        self._update_letter_type(letter_type)

    @classmethod
    def from_lettercode(cls, lettercode):
        """
        Initializes with the given LetterCode object.

        Parameters
        ----------
        lettercode : LetterCode
            LetterCode object.

        Returns
        -------
        LetterCode
            Copy of lettercode LetterCode object

        Raises
        ------
        TypeError
            If lettercode is of the wrong type.
        """
        if isinstance(lettercode, LetterCode):
            return cls(lettercode.letter_code, lettercode.letter_type)
        raise TypeError("lettercode must be a LetterCode")

    @property
    def letter_code(self):
        """return letter_code."""
        return self._letter_code

    @property
    def letter_type(self):
        """return letter_type."""
        return self._letter_type

    @letter_type.setter
    def letter_type(self, letter_type_value):
        """
        Sets letter_type and updates all other relevant properties as needed.

        Parameters
        ----------
        letter_type_value : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.

        Raises
        ------
        TypeError
            If letter_type_value is of the wrong type.
        """
        self._update_letter_type(letter_type_value)

    @letter_type.deleter
    def letter_type(self):
        """
        Sets letter_type to the default value (None) and updates all other relevant properties as needed.
        """
        self._update_letter_type(None)

    @property
    def description(self):
        """return description."""
        description = ""
        if self._letter_type in LETTER_CODES:
            if self._letter_code in LETTER_CODES[self._letter_type][0]:  # good
                description = LETTER_CODES[self._letter_type][0][self._letter_code]
            elif self._letter_code in LETTER_CODES[self._letter_type][1]:  # degenerate
                description = LETTER_CODES[self._letter_type][1][self._letter_code]
        return description

    @property
    def degenerate(self):
        """return degenerate."""
        return self._degenerate

    @property
    def supported(self):
        """return supported."""
        return self._supported

    @property
    def in_fasta_spec(self):
        """return in_fasta_spec."""
        return self._in_fasta_spec

    def complement(self):
        """
        Complementary letter code (ideally, of a nucleotide).

        If letter_code is not a nucleotide letter code, the complementary will be letter_code.
        In order not to impose the setting of letter_type as 'nucleotide', this method will work for any letter code
        (as long as letter_type is not 'aminoacid'), which has the side effect of returning nonsensical results when
        letter_code is not a nucleotide.
        Ex: For aminoacid letter codes that overlap with nucleotide letter codes, the output will be the complement of
        the nucleotide represented by the same letter code, which makes no sense.

        Returns
        -------
        LetterCode
            Complement of current LetterCode. Same LetterCode is returned if letter code is not a valid nucleotide.

        Raises
        ------
        TypeError
            If self.letter_type is 'aminoacid'.
        """
        if self._letter_type == "aminoacid":
            raise TypeError(
                "Complement is not possible for aminoacids (letter_type == 'aminoacid')"
            )
        if self._letter_type is None:
            warnings.warn(
                "letter_type is not explicitly 'nucleotide'. "
                "Therefore, the complementary letter code might not make sense."
            )
        return LetterCode(
            NUCLEOTIDE_LETTER_CODES_COMPLEMENT.get(
                self._letter_code, self._letter_code
            ),
            self._letter_type,
        )

    def _update_letter_type(self, letter_type):
        """
        Updates letter_type and all other relevant properties as needed.

        Parameters
        ----------
        letter_type : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.

        Raises
        ------
        TypeError
            If letter_type is of the wrong type.
        """
        if letter_type in LETTER_CODES:
            self._letter_type = letter_type

            if (
                self._letter_code in LETTER_CODES[self._letter_type][0]
            ):  # letter_codes_good
                self._degenerate = False
                self._supported = True
            elif (
                self._letter_code in LETTER_CODES[self._letter_type][1]
            ):  # letter_codes_degenerate
                self._degenerate = True
                self._supported = True
            else:  # _letter_code isn't defined in the FASTA specification
                self._degenerate = None
                self._supported = False
        elif letter_type is None:
            self._letter_type = letter_type
            self._supported = False
            self._degenerate = None
        else:
            raise TypeError(
                "letter_type must be one of: %s or None" % ", ".join(LETTER_CODES)
            )

    def __eq__(self, other):
        """
        Two LetterCode objects are equal if they represent the same letter code.
        A LetterCode is equal to a string if that string is the same as its letter code.
        """
        if isinstance(other, LetterCode):
            return self._letter_code == other.letter_code
        if isinstance(other, str):
            return self._letter_code == other.upper()
        return False

    def __repr__(self):
        return "LetterCode(%r)" % self._letter_code

    def __str__(self):
        return self._letter_code
