#!python
# coding: utf-8

"""
FastaSequence - Represents a single DNA/RNA/aminoacid FASTA sequence.
"""

import warnings
from .constants import LETTER_CODES, AMINOACIDS_NOT_IN_NUCLEOTIDES
from .lettercode import LetterCode


warnings.simplefilter(
    "always"
)  # show warnings everytime instead of only the first time they happen


class FastaSequence:
    """
    Represents one FASTA sequence.
    The class itself is an iterator of the sequence it represents.

    Attributes
    ----------
    id : str
        ID portion of the definition line (header). Can be empty.
    description : str
        Description portion of the definition line (header). Can be empty.
    sequence : list of LetterCode
        Sequence.
    sequence_type : 'nucleotide', 'aminoacid' or None
        Indicates the type of sequence ('aminoacid' or 'nucleotide'). Can be None if not known.
    inferred_type: bool
        True if FastaSequence inferred the sequence type, False otherwise.

    Methods
    -------
    from_fastasequence(fastasequence)
        Alternate __init__ method. Initializes instance with a FastaSequence object as only parameter.
    complement(reverse=False)
        Returns the complementary FastaSequence of a nucleotide sequence.
    gc_content(as_percentage=False)
        Returns the GC content of a nucleotide sequence.
    at_gc_ratio()
        Returns the AT/GC ratio of a nucleotide sequence.
    count_letter_codes(letter_codes=None)
        Returns the counts of each letter code if specified or all of the letter codes in the sequence if not.
    count_letter_codes_degenerate()
        Returns the counts of each degenerate letter code in the sequence.
    formatted_definition_line()
        Returns a formatted FASTA definition line (header).
    formatted_sequence(max_characters_per_line=70)
        Returns a formatted FASTA sequence (only the sequence, without the definition line).
    formatted_fasta()
        Returns a formatted FASTA (definition line and sequence).
    sequence_as_string()
        Returns the sequence as string.
    reverse()
        Iterator over the sequence of LetterCodes in reverse order.

    Raises
    ------
    TypeError
        When calling __init__, if sequence, id_, description, sequence_type or infer_type are of the wrong type.
        When calling from_fastasequence(), if fastasequence is of the wrong type.
        When setting id, if id_value is not str.
        When setting description, if description_value is not str.
        When setting sequence_type, if sequence_type_value is of the wrong type.
        When calling complement(), if sequence_type is 'aminoacid' or reverse is not bool.
        When calling gc_content(), if sequence_type is 'aminoacid' or as_percentage is not bool.
        When calling at_gc_ratio(), if sequence_type is 'aminoacid'.
        When calling count_letter_codes(), if letter_codes is not an iterable or None.
        When calling count_letter_codes_degenerate(), if self._sequence_type is not explicitly defined.
        When calling formatted_sequence(), if max_characters_per_line is not an int.
        When calling __getitem__, if item is not an int/slice or the sliced sequence is empty.
    """

    def __init__(
        self, sequence, id_="", description="", sequence_type=None, infer_type=False
    ):
        """
        Initializes FASTA sequence.

        Parameters
        ----------
        sequence : str
            String of characters representing a DNA, RNA or aminoacid sequence.
            Must be provided and cannot be empty.
        id_ : str, optional
            ID portion of the definition line (header).
            '>' and newlines will be removed, if any. Spaces will be converted to '_'. Can be an empty string.
        description : str, optional
            Description portion of the definition line (header).
            Newlines will be removed, if any. Can be an empty string.
        sequence_type : 'nucleotide', 'aminoacid' or None, optional
            Indicates the type of sequence ('aminoacid' or 'nucleotide').
            If not defined, FastaSequence can try to infer type based on the letter codes.
        infer_type : bool, optional
            Indicates if FastaSequence should try to infer aminoacid sequence type.
            If True, FastaSequence will analyse the whole sequence, in the worst case scenario,
            and can only identify aminoacid sequences.

        Raises
        ------
        TypeError
            If sequence, id_, description, sequence_type or infer_type are of the wrong type.
        """
        self._update_id(id_)
        self._update_description(description)
        self._update_sequence_type(sequence_type, update_letter_code_objects=False)

        if isinstance(sequence, str) and len(sequence) > 0:
            if isinstance(infer_type, bool):
                if infer_type:
                    self._sequence_type = self._infer_sequence_type(sequence)
                    # if infer_type is False there is no need to set _inferred_type as False
                    # as it is already set as such in _update_sequence_type
            else:
                raise TypeError("infer_type must be bool")

            # _sequence = [LetterCode, ...]
            # _counts = {letter: count, ...}
            self._sequence, self._counts = self._build_letter_code_sequence_and_counts(
                sequence
            )
        else:
            raise TypeError("sequence must be a non empty str")

        self._current_iterator = None
        self._gc = None
        self._at = None

    @classmethod
    def from_fastasequence(cls, fastasequence):
        """
        Initializes with the given FastaSequence object.

        Parameters
        ----------
        fastasequence : FastaSequence
            FastaSequence object.

        Returns
        -------
        FastaSequence
            Copy of fastasequence FastaSequence object.

        Raises
        ------
        TypeError
            If fastasequence is of the wrong type.
        """
        if isinstance(fastasequence, FastaSequence):
            return cls(
                fastasequence.sequence_as_string(),
                fastasequence.id,
                fastasequence.description,
                fastasequence.sequence_type,
            )
        raise TypeError("fastasequence must be a FastaSequence")

    @property
    def id(self):
        """return id."""
        return self._id

    @id.setter
    def id(self, id_value):
        """
        Sets id.
        '>' and newlines will be removed, if any. Spaces will be converted to '_'.
        Can be an empty string.

        Parameters
        ----------
        id_value : str
            ID portion of the definition line (header).

        Raises
        ------
        TypeError
            If id_value is not str.
        """
        self._update_id(id_value)

    @id.deleter
    def id(self):
        """
        Sets id to the default value ('').
        """
        self._update_id("")

    @property
    def description(self):
        """return description."""
        return self._description

    @description.setter
    def description(self, description_value):
        """
        Sets description.
        Newlines will be removed, if any.
        Can be an empty string.

        Parameters
        ----------
        description_value : str
            Description portion of the definition line (header).

        Raises
        ------
        TypeError
            If description_value is not str.
        """
        self._update_description(description_value)

    @description.deleter
    def description(self):
        """
        Sets description to the default value ('').
        """
        self._update_description("")

    @property
    def sequence(self):
        """return sequence."""
        return self._sequence

    @property
    def sequence_type(self):
        """return sequence_type."""
        return self._sequence_type

    @sequence_type.setter
    def sequence_type(self, sequence_type_value):
        """
        Sets sequence_type and updates all other relevant properties as needed.

        Parameters
        ----------
        sequence_type_value : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.

        Raises
        ------
        TypeError
            If sequence_type_value is of the wrong type.
        """
        self._update_sequence_type(sequence_type_value)

    @sequence_type.deleter
    def sequence_type(self):
        """
        Sets sequence_type to the default value (None) and updates all other relevant properties as needed.
        """
        self._update_sequence_type(None)

    @property
    def inferred_type(self):
        """return inferred_type."""
        return self._inferred_type

    def complement(self, reverse=False):
        """
        Complementary sequence (ideally, of a nucleotide sequence).

        Non-nucleotide LetterCodes don't have a complement and, therefore, stay the same.
        In order not to impose the setting of sequence_type as 'nucleotide', this method will work for any sequence and
        LetterCode (as long as sequence_type is not 'aminoacid'), which has the side effect of returning nonsensical
        results when LetterCodes are not nucleotides.
        Ex: For aminoacid letter codes that overlap with nucleotide letter codes, the output will be the complement of
        the nucleotide represented by the same letter code, which makes no sense.

        Parameters
        ----------
        reverse : bool, optional
            If sequence should be reversed.

        Returns
        -------
        FastaSequence
            Complement of the current nucleotide FastaSequence. Non-nucleotide LetterCodes will stay the same.

        Raises
        ------
        TypeError
            If self.sequence_type is 'aminoacid'.
            If reverse is not bool.
        """
        if isinstance(reverse, bool):
            if self._sequence_type == "aminoacid":
                raise TypeError(
                    "Complement is not possible for aminoacid sequences (sequence_type == 'aminoacid')"
                )
            if self._sequence_type is None:
                warnings.warn(
                    "sequence_type is not explicitly 'nucleotide'. "
                    "Therefore, the complementary sequence might not make sense."
                )
            if reverse:
                complement_sequence = "".join(
                    [
                        letter.complement().letter_code
                        for letter in reversed(self._sequence)
                    ]
                )
                reversed_text = "REVERSE "
            else:
                complement_sequence = "".join(
                    [letter.complement().letter_code for letter in self._sequence]
                )
                reversed_text = ""

            space = " " if len(self._description) > 0 else ""
            complement_description = "%s[%sCOMPLEMENT]" % (space, reversed_text)
            return FastaSequence(
                complement_sequence,
                self._id,
                self._description + complement_description,
                self._sequence_type,
            )
        raise TypeError("reverse must be a bool")

    def gc_content(self, as_percentage=False):
        """
        Calculates the GC content of nucleotide sequence (as a ratio, by default).
        Ignores degenerate letter codes besides S (G or C).
        GC content is calculated the first time the method is called. Later calls will retrieve the same value.
        GC content can also be calculated in at_gc_ratio.
        If sequence_type is not 'nucleotide' (or the sequence is not inherently a nucleotide sequence) the GC content
        might be nonsensical.

        Parameters
        ----------
        as_percentage : bool, optional
            Indicates whether the computed value should be returned as a percentage instead of the default ratio.

        Returns
        -------
        float
            GC content of sequence.

        Raises
        ------
        TypeError
            If self.sequence_type is 'aminoacid'.
            If as_percentage is not bool.
        """
        if self._sequence_type == "aminoacid":
            raise TypeError(
                "GC content is not meant to be calculated for aminoacid sequences "
                "(sequence_type == 'aminoacid')"
            )
        if self._sequence_type is None:
            warnings.warn(
                "sequence_type is not explicitly 'nucleotide'. "
                "Therefore, the calculated GC content might not make sense."
            )
        if isinstance(as_percentage, bool):
            if not self._gc:  # if gc_content was not called before
                gc = 0
                for letter_code in self._sequence:
                    if letter_code.letter_code in (
                        "G",
                        "C",
                        "S",
                    ):  # S means either G or C
                        gc += 1
                self._gc = gc
            gc_content = self._gc / len(self._sequence)
            return gc_content * 100 if as_percentage else gc_content
        raise TypeError("as_percentage must be a bool")

    def at_gc_ratio(self):
        """
        Calculates the AT/GC ratio of nucleotide sequence.
        Ignores degenerate letter codes besides W (A or T) and S (G or C).
        AT/GC ratio is calculated the first time the method is called. Later calls will retrieve the same value.
        Also uses previously calculated _gc or calculates and saves it if it hasn't been calculated yet.
        If sequence_type is not 'nucleotide' (or the sequence is not inherently a nucleotide sequence) the AT/GC ratio
        might be nonsensical.

        Returns
        -------
        float
            AT/GC ratio of sequence.

        Raises
        ------
        TypeError
            If self.sequence_type is 'aminoacid'.
        """
        if self._sequence_type == "aminoacid":
            raise TypeError(
                "AT/GC ratio is not meant to be calculated for aminoacid sequences "
                "(sequence_type == 'aminoacid')"
            )
        if self._sequence_type is None:
            warnings.warn(
                "sequence_type is not explicitly 'nucleotide'. "
                "Therefore, the calculated AT/GC ratio might not make sense."
            )
        if not self._at or not self._gc:  # if at_gc_ratio was not called before
            at = 0
            gc = 0
            for letter_code in self._sequence:
                if self._at is None and letter_code.letter_code in (
                    "A",
                    "T",
                    "W",
                ):  # W means either A or T
                    at += 1
                elif self._gc is None and letter_code.letter_code in (
                    "G",
                    "C",
                    "S",
                ):  # S means either G or C
                    gc += 1
            if self._gc is None:
                self._gc = gc
            if self._at is None:
                self._at = at
        return self._at / self._gc if self._gc != 0 else 0

    def count_letter_codes(self, letter_codes=None):
        """
        Returns letter code counts.
        By default shows counts for all existing letter codes in the sequence,
        but specific letter codes can be specified.

        Parameters
        ----------
        letter_codes : iterable or None, optional
            Iterable of all letter codes to count.

        Returns
        -------
        dict
            Counts for every letter code in letter_codes or all letter codes in the sequence
            if letter_codes is not specified.

        Raises
        ------
        TypeError
            If letter_codes is not an iterable or None.
        """
        if letter_codes is None or not letter_codes:
            return self._counts
        return {letter: self._counts.get(letter, 0) for letter in iter(letter_codes)}

    def count_letter_codes_degenerate(self):
        """
        Returns degenerate letter code counts.
        sequence_type must be explicitly defined for this method to work.

        Returns
        -------
        dict
            Counts for every degenerate letter code in the sequence.

        Raises
        ------
        TypeError
            If self._sequence_type is not explicitly defined.
        """
        if self._sequence_type in LETTER_CODES:
            return {
                letter: counts
                for letter, counts in self._counts.items()
                if letter in LETTER_CODES[self._sequence_type][1]
            }

        raise TypeError(
            "To count degenerate letter codes the sequence_type must be explicitly '%s'"
            % "' or '".join(LETTER_CODES)
        )

    def formatted_definition_line(self):
        """
        Formatted FASTA definition line (header).

        Returns
        -------
        str
            FASTA definition line properly formatted.
        """
        return (
            ">%s %s" % (self._id, self._description)
            if self._description
            else ">" + self._id
        )

    def formatted_sequence(self, max_characters_per_line=70):
        """
        Formatted FASTA sequence (only the sequence, without the definition line).
        Lines are separated by '\n'.

        Parameters
        ----------
        max_characters_per_line : int
            Maximum number of characters per line.
            This value should not go above 80, as per the FASTA specification.
            A very low value is also not recommended.

        Returns
        -------
        str
            FASTA sequence properly formatted.

        Raises
        ------
        TypeError
            If max_characters_per_line is not an int.
        """
        if isinstance(max_characters_per_line, int):
            max_characters_per_line = (
                1 if max_characters_per_line <= 0 else max_characters_per_line
            )

            current_character_count = 0
            final_sequence = ""
            while current_character_count < len(self._sequence):
                temp_sequence = self._sequence[
                    current_character_count : current_character_count
                    + max_characters_per_line
                ]
                final_sequence += "".join(map(str, temp_sequence)) + "\n"
                current_character_count += max_characters_per_line
            return final_sequence[:-1]  # remove last '\n'
        raise TypeError("max_characters_per_line must be an int")

    def formatted_fasta(self):
        """
        Formatted FASTA (definition line and sequence).

        Returns
        -------
        str
            FASTA properly formatted.
        """
        return self.formatted_definition_line() + "\n" + self.formatted_sequence()

    def sequence_as_string(self):
        """
        Returns the sequence as string.
        Converts the list of LetterCode objects to a single string.

        Returns
        -------
        str
            Sequence as string.
        """
        return "".join(map(str, self._sequence))

    def reverse(self):
        """
        Iterates over the sequence in reverse order (same as calling reversed() on a FastaSequence object).
        Returns a new iterator of the sequence (from the end) every time reverse is called.

        Returns
        -------
        iterator
            Iterator over the reversed sequence.
        """
        return reversed(self)

    def _update_id(self, id_):
        """
        Updates ID portion of the definition line (header).
        '>' and newlines will be removed, if any. Spaces will be converted to '_'.
        Can be an empty string.

        Parameters
        ----------
        id_ : str
            ID portion of the definition line (header).

        Raises
        ------
        TypeError
            If id_ is not str.
        """
        if isinstance(id_, str):
            id_ = "".join(
                id_.strip().replace(" ", "_").split()
            )  # remove spaces and newlines
            if id_.startswith(">"):  # remove '>' if any
                id_ = id_[1:]
            self._id = id_
        else:
            raise TypeError("id_ must be str")

    def _update_description(self, description):
        """
        Updates description portion of the definition line (header).
        Newlines will be removed, if any.
        Can be an empty string.

        Parameters
        ----------
        description : str
            Description portion of the definition line (header).

        Raises
        ------
        TypeError
            If description is not str.
        """
        if isinstance(description, str):
            self._description = " ".join(
                description.strip().split()
            )  # remove extra spaces and newlines
        else:
            raise TypeError("description must be str")

    def _update_sequence_type(self, sequence_type, update_letter_code_objects=True):
        """
        Updates sequence_type and all other relevant properties as needed.

        Parameters
        ----------
        sequence_type : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.
        update_letter_code_objects : bool
            Should LetterCode objects be updated with sequence_type or not.

        Raises
        ------
        TypeError
            If sequence_type or update_letter_code_objects are of the wrong type.
        """
        if (
            isinstance(sequence_type, str) and sequence_type in LETTER_CODES
        ) or sequence_type is None:
            self._sequence_type = sequence_type
            self._inferred_type = False
        else:
            raise TypeError(
                "sequence_type must be one of: '%s' or None" % "', '".join(LETTER_CODES)
            )
        if isinstance(update_letter_code_objects, bool):
            if update_letter_code_objects:
                for letter_code_object in self._sequence:  # update LetterCode objects
                    letter_code_object.letter_type = self._sequence_type
        else:
            raise TypeError("update_letter_code_objects must be a bool")

    def _build_letter_code_sequence_and_counts(self, string_sequence):
        """
        Iterates over the sequence and builds a list of LetterCode objects while counting the number of letter codes.

        Parameters
        ----------
        string_sequence: str
            String of characters representing a DNA, RNA or aminoacid sequence.

        Returns
        -------
        (list of LetterCode, dict of letter code counts)
        """
        letter_code_list = []
        letter_code_count_dict = {}
        for letter_code in string_sequence:
            letter_code = letter_code.upper()
            letter_code_list.append(LetterCode(letter_code, self._sequence_type))
            letter_code_count_dict[letter_code] = (
                letter_code_count_dict.setdefault(letter_code, 0) + 1
            )

        return letter_code_list, letter_code_count_dict

    def _infer_sequence_type(self, string_sequence):
        """
        Tries to infer aminoacid sequence type.
        Tests for the presence of letter codes that can only represent aminoacids
        (ie, aminoacids letter codes not in nucleotides letter codes).
        The reverse (testing for nucleotides) is not 100% accurate because there are no letter codes
        which belong solely to nucleotide type sequences.
        Assumes no unknown letter codes.

        Parameters
        ----------
        string_sequence: str
            String of characters representing a DNA, RNA or aminoacid sequence.

        Returns
        -------
        'aminoacid' or existing value (can be None)
            Inferred type 'aminoacid', existing value otherwise.
        """
        for letter_code in string_sequence:
            if letter_code in AMINOACIDS_NOT_IN_NUCLEOTIDES:
                self._inferred_type = True
                return "aminoacid"
        # self._inferred_type = False  # _inferred_type is already False when this function is called in __init__
        return self._sequence_type  # returns the already set value

    def __iter__(self):
        """
        Iterates over the sequence.
        Returns a new iterator of the sequence (from the beginning) every time __iter__ is called.
        """

        def iter_sequence():
            for letter_code in self._sequence:
                yield letter_code

        self._current_iterator = iter_sequence()
        return self._current_iterator

    def __reversed__(self):
        """
        Iterates over the sequence in reverse.
        Returns a new iterator of the reversed sequence (from the end) every time __reversed__ is called.
        """

        def iter_sequence_reversed():
            for letter_code in reversed(self._sequence):
                yield letter_code

        self._current_iterator = iter_sequence_reversed()
        return self._current_iterator

    def __next__(self):
        """
        Returns the next letter code from the current iterator (most recent iterator).
        If no iterator still exists, calls __iter__ to create it.
        """
        if self._current_iterator is None:
            self.__iter__()
        return next(self._current_iterator)

    def __getitem__(self, item):
        """
        Indexing returns the LetterCode object at the given index.
        Slicing returns a new FastaSequence copy of self with the sliced sequence of LetterCode objects.
        The description line is updated to reflect the slices made to the original sequence.
        Slices can't return an empty sequence.

        Parameters
        ----------
        item : int or slice

        Returns
        -------
        LetterCode, FastaSequence
            LetterCode at given index.
            or
            FastaSequence with a sliced sequence of LetterCode objects.
        """
        if isinstance(item, int):
            return self._sequence[item]
        if isinstance(item, slice):
            new_sequence = self.sequence_as_string()[item]
            if len(new_sequence) == 0:
                raise TypeError(
                    "Slice resulted in an empty sequence. FastaSequence must have a non-empty sequence"
                )
            slice_text = "[SLICE OF ORIGINAL: %s:%s:%s]" % (
                item.start,
                item.stop,
                item.step,
            )
            new_description = (
                "%s %s" % (self.description, slice_text)
                if self.description
                else slice_text
            )
            return FastaSequence(
                new_sequence, self.id, new_description, self.sequence_type
            )
        raise TypeError("Indices must be integers or slices")

    def __eq__(self, other):
        """
        Two FastaSequence objects are equal if they represent the same sequence.
        A FastaSequence is equal to a string if it represents the same string sequence.
        A FastaSequence is equal to a list if it represents the same LetterCode sequence.
        """
        if isinstance(other, FastaSequence):
            return self._sequence == other.sequence
        if isinstance(other, str):
            return self.sequence_as_string() == other
        if isinstance(other, list):
            return self._sequence == other
        return False

    def __len__(self):
        return len(self._sequence)

    def __repr__(self):
        return "FastaSequence(%r)" % self.sequence_as_string()

    def __str__(self):
        """
        >id description
        sequence
        """
        return self.formatted_fasta()
