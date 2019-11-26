#!python
# coding: utf-8

"""
FastaParser

Parses FASTA files with the Reader class (generates FastaSequence objects or plane strings of the parsed sequences)
Writes FASTA files with the Writer class (takes FastaSequence objects or headers + sequences as strings)

ex:
    > import FastaParser
    > with open('fasta_file.fasta') as fasta_file:
    >   reader = FastaParser.Reader(fasta_file)
    >   [seq.id for seq in reader]
    ['HSBGPG', 'HSGLTH1']

ex:
    > import FastaParser
    > with open('fasta_file.fasta', 'w') as fasta_file:
    >   writer = FastaParser.Writer(fasta_file)
    >   seqs = [('HSBGPG example sequence', 'TTCCAGGTGTGCCAATCCAGTCCATG'),
    >   ...     ('HSGLTH1 example sequence 2', 'GTACCTGACCTAACCGTGTGGACCTT')]
    > writer.writefastas(seqs)

Based on these pages:
http://genetics.bwh.harvard.edu/pph/FASTA.html
https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
https://en.wikipedia.org/wiki/FASTA_format
https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta
http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
"""


import warnings
from collections import namedtuple


__author__ = 'Pedro HC David, https://github.com/Kronopt'
__credits__ = ['Pedro HC David']
__version__ = '0.1'


warnings.simplefilter("always")  # show warnings everytime instead of only the first time they happen

# NUCLEOTIDE DICTIONARIES
NUCLEOTIDE_LETTER_CODES_GOOD = {
    'A': 'adenosine',
    'C': 'cytidine',
    'G': 'guanine',
    'T': 'thymidine',
    'N': 'any (A/G/C/T)',
    'U': 'uridine'
}
NUCLEOTIDE_LETTER_CODES_DEGENERATE = {
    'K': 'keto (G/T)',
    'S': 'strong (G/C)',
    'Y': 'pyrimidine (T/C)',
    'M': 'amino (A/C)',
    'W': 'weak (A/T)',
    'R': 'purine (G/A)',
    'B': 'G/T/C',
    'D': 'G/A/T',
    'H': 'A/C/T',
    'V': 'G/C/A',
    '-': 'gap of indeterminate length'
}
NUCLEOTIDE_LETTER_CODES_COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
    'U': 'A',
    'K': 'M',
    'S': 'S',
    'Y': 'R',
    'M': 'K',
    'W': 'W',
    'R': 'Y',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
    '-': '-'
}

# AMINOACID DICTIONARIES
AMINOACID_LETTER_CODES_GOOD = {
    'A': 'alanine',
    'B': 'aspartate/asparagine',
    'C': 'cystine',
    'D': 'aspartate',
    'E': 'glutamate',
    'F': 'phenylalanine',
    'G': 'glycine',
    'H': 'histidine',
    'I': 'isoleucine',
    'K': 'lysine',
    'L': 'leucine',
    'M': 'methionine',
    'N': 'asparagine',
    'P': 'proline',
    'Q': 'glutamine',
    'R': 'arginine',
    'S': 'serine',
    'T': 'threonine',
    'U': 'selenocysteine',
    'V': 'valine',
    'W': 'tryptophan',
    'Y': 'tyrosine',
    'Z': 'glutamate/glutamine',
    'X': 'any',
    '*': 'translation stop'
}
AMINOACID_LETTER_CODES_DEGENERATE = {
    '-': 'gap of indeterminate length'
}

# BOTH
LETTER_CODES = {
    'nucleotide': (NUCLEOTIDE_LETTER_CODES_GOOD, NUCLEOTIDE_LETTER_CODES_DEGENERATE),
    'aminoacid': (AMINOACID_LETTER_CODES_GOOD, AMINOACID_LETTER_CODES_DEGENERATE)
}

# set operations
LETTER_CODES_ALL = set(
    list(NUCLEOTIDE_LETTER_CODES_GOOD) +
    list(NUCLEOTIDE_LETTER_CODES_DEGENERATE) +
    list(AMINOACID_LETTER_CODES_GOOD) +
    list(AMINOACID_LETTER_CODES_DEGENERATE)
)
NUCLEOTIDE_LETTER_CODES_ALL = set(
    list(NUCLEOTIDE_LETTER_CODES_GOOD) +
    list(NUCLEOTIDE_LETTER_CODES_DEGENERATE)
)
AMINOACID_LETTER_CODES_ALL = set(
    list(AMINOACID_LETTER_CODES_GOOD) +
    list(AMINOACID_LETTER_CODES_DEGENERATE)
)
AMINOACIDS_NOT_IN_NUCLEOTIDES = AMINOACID_LETTER_CODES_ALL - NUCLEOTIDE_LETTER_CODES_ALL
# nucleotides_not_in_aminoacids would be empty


class LetterCode:
    """
    Represents a single letter code.

    Attributes
    ----------
    letter_code : str
        Upper case letter code.
    sequence_type : str or None
        'nucleotide' or 'aminoacid. None if there is no information about sequence type.
    description : str
        Description or nucleotide/aminoacid name of letter code (can be an empty string).
    degenerate : bool or None
        Indicates if a letter code is degenerate or not (can be None if letter code is
        not defined in the FASTA specification or sequence_type is unknown).
    supported : bool
        Indicates if letter code is supported or not
        (ie, if sequence_type is provided and letter code is defined in the FASTA specification).
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
        When calling __init__, if letter_code or sequence_type are of the wrong type.
        When calling from_lettercode(), if lettercode is of the wrong type.
        When setting sequence_type, if sequence_type_value is of the wrong type.
        When calling complement(), if sequence_type is 'aminoacid'.
    """

    def __init__(self, letter_code, sequence_type=None):
        """
        Initializes letter code given.

        Parameters
        ----------
        letter_code : str
            Letter code.
        sequence_type : 'nucleotide', 'aminoacid' or None, optional
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.

        Raises
        ------
        TypeError
            If letter_code or sequence_type are of the wrong type.
        """
        if isinstance(letter_code, str) and len(letter_code) == 1:
            self._letter_code = letter_code.upper()
            self._in_fasta_spec = self._letter_code in LETTER_CODES_ALL
        else:
            raise TypeError('letter_code must be a single character str')

        self._update_sequence_type(sequence_type)

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
            return cls(lettercode.letter_code, lettercode.sequence_type)
        else:
            raise TypeError('lettercode must be a LetterCode')

    @property
    def letter_code(self):
        return self._letter_code

    @property
    def sequence_type(self):
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

    @property
    def description(self):
        description = ''
        if self._sequence_type in LETTER_CODES:
            if self._letter_code in LETTER_CODES[self._sequence_type][0]:  # good
                description = LETTER_CODES[self._sequence_type][0][self._letter_code]
            elif self._letter_code in LETTER_CODES[self._sequence_type][1]:  # degenerate
                description = LETTER_CODES[self._sequence_type][1][self._letter_code]
        return description

    @property
    def degenerate(self):
        return self._degenerate

    @property
    def supported(self):
        return self._supported

    @property
    def in_fasta_spec(self):
        return self._in_fasta_spec

    def complement(self):
        """
        Complementary letter code (ideally, of a nucleotide).

        If letter_code is not a nucleotide letter code, the complementary will be letter_code.
        In order not to impose the setting of sequence_type as 'nucleotide', this method will work for any letter code
        (as long as sequence_type is not 'aminoacid'), which has the side effect of returning nonsensical results when
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
            If self.sequence_type is 'aminoacid'.
        """
        if self._sequence_type == 'aminoacid':
            raise TypeError('Complement is not possible for aminoacids (sequence_type == \'aminoacid\')')
        if self._sequence_type is None:
            warnings.warn('sequence_type is not explicitly \'nucleotide\'. '
                          'Therefore, the complementary letter code might not make sense.')
        return LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT.get(self._letter_code, self._letter_code),
                          self._sequence_type)

    def _update_sequence_type(self, sequence_type):
        """
        Updates sequence_type and all other relevant properties as needed.

        Parameters
        ----------
        sequence_type : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.

        Raises
        ------
        TypeError
            If sequence_type is of the wrong type.
        """
        if sequence_type in LETTER_CODES:
            self._sequence_type = sequence_type

            if self._letter_code in LETTER_CODES[self._sequence_type][0]:  # letter_codes_good
                self._degenerate = False
                self._supported = True
            elif self._letter_code in LETTER_CODES[self._sequence_type][1]:  # letter_codes_degenerate
                self._degenerate = True
                self._supported = True
            else:  # _letter_code isn't defined in the FASTA specification
                self._degenerate = None
                self._supported = False
        elif sequence_type is None:
            self._sequence_type = sequence_type
            self._supported = False
            self._degenerate = None
        else:
            raise TypeError('sequence_type must be one of: %s or None' % ', '.join(LETTER_CODES))

    def __eq__(self, other):
        """
        Two LetterCode objects are equal if they represent the same letter code.
        A LetterCode is equal to a string if that string is the same as its letter code.
        """
        if isinstance(other, LetterCode):
            return self._letter_code == other.letter_code
        elif isinstance(other, str):
            return self._letter_code == other.upper()
        else:
            return False

    def __repr__(self):
        return 'LetterCode(%r)' % self._letter_code

    def __str__(self):
        return self._letter_code


# TODO Identify FASTA ID's (see linked sources)
# TODO FASTAID class (?) to then return in the id property
# TODO FASTQ parser
# TODO per fasta sequence, show warning if there are characters not in the FASTA specification
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
        When setting sequence_type, if sequence_type_value is of the wrong type.
        When calling complement(), if sequence_type is 'aminoacid' or reverse is not bool.
        When calling gc_content(), if sequence_type is 'aminoacid' or as_percentage is not bool.
        When calling at_gc_ratio(), if sequence_type is 'aminoacid'.
        When calling count_letter_codes(), if letter_codes is not an iterable or None.
        When calling count_letter_codes_degenerate(), if self._sequence_type is not explicitly defined.
        When calling formatted_sequence(), if max_characters_per_line is not an int.
        When calling __getitem__, if item is not an int/slice or the sliced sequence is empty.
    """

    def __init__(self, sequence, id_='', description='', sequence_type=None, infer_type=False):
        """
        Initializes FASTA sequence.

        Parameters
        ----------
        sequence : str
            String of characters representing a DNA, RNA or aminoacid sequence.
            Must be provided and cannot be empty.
        id_ : str, optional
            ID portion of the definition line (header).
            Can be an empty string.
        description : str, optional
            Description portion of the definition line (header).
            Can be an empty string.
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
        if isinstance(id_, str):
            self._id = id_
        else:
            raise TypeError('id_ must be str')

        if isinstance(description, str):
            self._description = description
        else:
            raise TypeError('description must be str')

        self._update_sequence_type(sequence_type, update_letter_code_objects=False)

        if isinstance(sequence, str) and len(sequence) > 0:
            if isinstance(infer_type, bool):
                if infer_type:
                    self._sequence_type = self._infer_sequence_type(sequence)
                    # if infer_type is False there is no need to set _inferred_type as False
                    # as it is already set as such in _update_sequence_type
            else:
                raise TypeError('infer_type must be bool')

            # _sequence = [LetterCode, ...]
            # _counts = {letter: count, ...}
            self._sequence, self._counts = self._build_letter_code_sequence_and_counts(sequence)
        else:
            raise TypeError('sequence must be a non empty str')

        self._gc_content = None
        self._at_gc_ratio = None

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
            Copy of fastasequence FastaSequence object

        Raises
        ------
        TypeError
            If fastasequence is of the wrong type.
        """
        if isinstance(fastasequence, FastaSequence):
            return cls(fastasequence.sequence_as_string(),
                       fastasequence.id,
                       fastasequence.description,
                       fastasequence.sequence_type)
        else:
            raise TypeError('fastasequence must be a FastaSequence')

    @property
    def id(self):
        return self._id

    @property
    def description(self):
        return self._description

    @property
    def sequence(self):
        return self._sequence

    @property
    def sequence_type(self):
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

    @property
    def inferred_type(self):
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
            if self._sequence_type == 'aminoacid':
                raise TypeError('Complement is not possible for aminoacid sequences (sequence_type == \'aminoacid\')')
            if self._sequence_type is None:
                warnings.warn('sequence_type is not explicitly \'nucleotide\'. '
                              'Therefore, the complementary sequence might not make sense.')
            if reverse:
                complement_sequence = ''.join([letter.complement().letter_code for letter in reversed(self._sequence)])
            else:
                complement_sequence = ''.join([letter.complement().letter_code for letter in self._sequence])

            return FastaSequence(complement_sequence, self._id, self._description, self._sequence_type)
        else:
            raise TypeError('reverse must be a bool')

    def gc_content(self, as_percentage=False):
        """
        Calculates the GC content of nucleotide sequence (as a ratio, by default).
        Ignores degenerate letter codes besides S (G or C).
        GC content is calculated the first time the method is called. Later calls will retrieve the same value.
        GC content can also be calculates in at_gc_ratio.
        If sequence_type is not 'nucleotide' (or the sequence is not inherently a nucleotide sequence) the GC content
        might be nonsensical.

        Parameters
        ----------
        as_percentage : bool, optional
            Indicates whether the computed value should be returned as a percentage instead of the default ratio.

        Returns
        -------
        float
            GC content of sequence

        Raises
        ------
        TypeError
            If self.sequence_type is 'aminoacid'.
            If as_percentage is not bool.
        """
        if isinstance(as_percentage, bool):
            if not self._gc_content:  # if gc_content was not called before
                if self._sequence_type == 'aminoacid':
                    raise TypeError('GC content is not meant to be calculated for aminoacid sequences '
                                    '(sequence_type == \'aminoacid\')')
                if self._sequence_type is None:
                    warnings.warn('sequence_type is not explicitly \'nucleotide\'. '
                                  'Therefore, the calculated GC content might not make sense.')
                gc = 0
                for letter_code in self._sequence:
                    if letter_code.letter_code in ('G', 'C', 'S'):  # S means either G or C
                        gc += 1
                self._gc_content = gc / (len(self._sequence))
            return self._gc_content * 100 if as_percentage else self._gc_content
        else:
            raise TypeError('as_percentage must be a bool')

    def at_gc_ratio(self):
        """
        Calculates the AT/GC ratio of nucleotide sequence.
        Ignores degenerate letter codes besides W (A or T) and S (G or C).
        AT/GC ratio is calculated the first time the method is called. Later calls will retrieve the same value.
        Also uses previously calculated _gc_content or calculates and saves it if it hasn't been calculated yet.
        If sequence_type is not 'nucleotide' (or the sequence is not inherently a nucleotide sequence) the AT/GC ratio
        might be nonsensical.

        Returns
        -------
        float
            AT/GC ratio of sequence

        Raises
        ------
        TypeError
            If self.sequence_type is 'aminoacid'.
        """
        if not self._at_gc_ratio:  # if at_gc_ratio was not called before
            if self._sequence_type == 'aminoacid':
                raise TypeError('AT/GC ratio is not meant to be calculated for aminoacid sequences '
                                '(sequence_type == \'aminoacid\')')
            if self._sequence_type is None:
                warnings.warn('sequence_type is not explicitly \'nucleotide\'. '
                              'Therefore, the calculated AT/GC ratio might not make sense.')
            at = 0
            gc = 0
            for letter_code in self._sequence:
                if letter_code.letter_code in ('A', 'T', 'W'):  # W means either A or T
                    at += 1
                elif self._gc_content is None and letter_code.letter_code in ('G', 'C', 'S'):  # S means either G or C
                    gc += 1
            if self._gc_content is None:
                self._gc_content = gc
            self._at_gc_ratio = at/self._gc_content if self._gc_content != 0 else 0
        return self._at_gc_ratio

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
            if letter_codes is not specified

        Raises
        ------
        TypeError
            If letter_codes is not an iterable or None.
        """
        if letter_codes is None or not letter_codes:
            return self._counts
        else:
            return {letter: self._counts.get(letter, 0) for letter in iter(letter_codes)}

    def count_letter_codes_degenerate(self):
        """
        Returns degenerate letter code counts.
        sequence_type must be explicitly defined for this method to work.

        Returns
        -------
        dict
            Counts for every degenerate letter code in the sequence

        Raises
        ------
        TypeError
            If self._sequence_type is not explicitly defined.
        """
        if self._sequence_type in LETTER_CODES:
            return {letter: counts for letter, counts in self._counts.items()
                    if letter in LETTER_CODES[self._sequence_type][1]}
        else:
            raise TypeError('To count degenerate letter codes the sequence_type must be '
                            'explicitly \'%s\'' % ('\' or \''.join(LETTER_CODES),))

    def formatted_definition_line(self):
        """
        Formatted FASTA definition line (header).

        Returns
        -------
        str
            FASTA definition line properly formatted
        """
        return '>' + self._id + ' ' + self._description

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
            FASTA sequence properly formatted

        Raises
        ------
        TypeError
            If max_characters_per_line is not an int.
        """
        if isinstance(max_characters_per_line, int):
            max_characters_per_line = 1 if max_characters_per_line <= 0 else max_characters_per_line

            current_character_count = 0
            final_sequence = ''
            while current_character_count < len(self._sequence):
                temp_sequence = self._sequence[current_character_count: current_character_count
                                               + max_characters_per_line]
                final_sequence += ''.join(map(str, temp_sequence)) + '\n'
                current_character_count += max_characters_per_line
            return final_sequence[:-1]  # remove last '\n'

        else:
            raise TypeError('max_characters_per_line must be an int')

    def formatted_fasta(self):
        """
        Formatted FASTA (definition line and sequence).

        Returns
        -------
        str
            FASTA properly formatted
        """
        return self.formatted_definition_line() + '\n' + self.formatted_sequence()

    def sequence_as_string(self):
        """
        Returns the sequence as string.
        Converts the list of LetterCode objects to a single string.

        Returns
        -------
        str
            Sequence as string
        """
        return ''.join(map(str, self._sequence))

    def reverse(self):
        """
        Iterates over the sequence in reverse order.
        Returns a new iterator of the sequence (from the end) every time reverse is called.

        Returns
        -------
        iterator
            Iterator over the reversed sequence
        """
        return self.__reversed__()

    def _update_sequence_type(self, sequence_type, update_letter_code_objects=True):
        """
        Updates sequence_type and all other relevant properties as needed.

        Parameters
        ----------
        sequence_type : 'nucleotide', 'aminoacid' or None
            'nucleotide' or 'aminoacid' type sequence, None if there is no information.
        update_letter_code_objects : bool
            Should LetterCode objects be updated with sequence_type or not

        Raises
        ------
        TypeError
            If sequence_type or update_letter_code_objects are of the wrong type.
        """
        if (isinstance(sequence_type, str) and sequence_type in LETTER_CODES) or sequence_type is None:
            self._sequence_type = sequence_type
            self._inferred_type = False
        else:
            raise TypeError('sequence_type must be one of: %s or None' % LETTER_CODES)
        if isinstance(update_letter_code_objects, bool):
            if update_letter_code_objects:
                for letter_code_object in self._sequence:  # update LetterCode objects
                    letter_code_object.sequence_type = self._sequence_type
        else:
            raise TypeError('update_letter_code_objects must be a bool')

    def _build_letter_code_sequence_and_counts(self, string_sequence):
        """
        Iterate over the sequence and build a list of LetterCode objects while counting the number of letter codes.

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
            letter_code_count_dict[letter_code] = letter_code_count_dict.setdefault(letter_code, 0) + 1

        return letter_code_list, letter_code_count_dict

    def _infer_sequence_type(self, string_sequence):
        """
        Tries to infer aminoacid sequence type.
        Tests for the presence of letter codes that can only represent aminoacids
        (ie, aminoacids letter codes not in nucleotides letter codes).
        The reverse (testing for nucleotides) is not 100% accurate because there are no letter codes
        which belong solely to nucleotide type sequences.

        Parameters
        ----------
        string_sequence: str
            String of characters representing a DNA, RNA or aminoacid sequence.

        Returns
        -------
        'aminoacid' or existing value (can be None)
            Inferred type 'aminoacid', existing value otherwise
        """
        for letter_code in string_sequence:
            if letter_code in AMINOACIDS_NOT_IN_NUCLEOTIDES:
                self._inferred_type = True
                return 'aminoacid'
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
        if not hasattr(self, '_current_iterator'):
            self.__iter__()
        return next(self._current_iterator)

    def __getitem__(self, item):
        """
        Slicing/indexing returns a new FastaSequence copy of self with the sliced sequence of LetterCode objects.
        The description line is updated to reflect the slices made to the original sequence.
        Slices can't return an empty sequence.

        Parameters
        ----------
        item : int or slice

        Returns
        -------
        FastaSequence
            Copy of self with the sliced sequence of LetterCode objects
        """
        if isinstance(item, (int, slice)):
            new_sequence = self.sequence_as_string()[item]
            if len(new_sequence) == 0:
                raise TypeError('Slice resulted in an empty sequence. FastaSequence must have a non-empty sequence')

            new_description = '%s [SLICE OF ORIGINAL: %s]' % (self.description, item)
            return FastaSequence(new_sequence, self.id, new_description, self.sequence_type)
        else:
            raise TypeError('Indices must be integers or slices')

    def __len__(self):
        return len(self._sequence)

    def __repr__(self):
        return 'FastaSequence(%r)' % self.sequence_as_string()

    def __str__(self):
        """
        >id description
        sequence
        """
        return ">%s %s\n%s" % (self._id, self._description, self.sequence_as_string())


###################
# READER AND WRITER
###################


class ParseDefinitionLine:
    """
    Implements a parser of FASTA definition lines (to be used by the Reader and Writer classes)

    Methods
    -------
    _parse_definition_line(definition_line)
        Parses FASTA definition lines

    Raises
    ------
    TypeError
        When calling _parse_definition_line, if definition_line is of the wrong type.
    """

    def _parse_definition_line(self, definition_line):
        """
        Parses a FASTA definition line and returns an id and a description.

        Parameters
        ----------
        definiton_line : str
            FASTA sequence definition line (header). May contain or not the '>' symbol at the start.
            Can be an empty string.

        Returns
        -------
        tuple
            ID and description. Can both be empty strings.

        Raises
        ------
        TypeError
            If definition_line is of the wrong type.
        """
        if isinstance(definition_line, str):  # '>id|more_id description ...' with or without the '>' at the start
            id_and_description = definition_line.split(maxsplit=1)  # first space separates id from description

            # both id and description can be empty
            if len(id_and_description) == 0 or (len(id_and_description) == 1 and id_and_description[0] == '>'):
                _id = ''
                _description = ''
            else:
                if id_and_description[0].startswith('>'):
                    id_and_description[0] = id_and_description[0][1:]
                if len(id_and_description) == 1:  # description can be empty (assumes only id present if len == 1)
                    _id = id_and_description[0]
                    _description = ''
                else:
                    _id, _description = id_and_description
        else:
            raise TypeError('definition_line must be str')

        return _id, _description


class Reader(ParseDefinitionLine):
    """
    Parser/Reader for the given FASTA file.
    Iterates over the FASTA file using one of two parsing mechanisms:
        'quick':
            Generates objects containing just the FASTA header and sequence attributes
            for each sequence in the FASTA file.
            Parses the FASTA file faster but lacks some features.
        'rich':
            Returns FastaSequence objects (default).
            Slower, but feature rich.

    Attributes
    ----------
    fasta_file : file object
        The FASTA file passed as parameter.
    sequences_type : 'nucleotide', 'aminoacid' or None
        Indicates the type of sequences to expect ('aminoacid' or 'nucleotide'). Can be None if not known.
    parse_method: 'rich' or 'quick'
        Parse method used ('rich' or 'quick').

    Raises
    ------
    TypeError
        When calling __init__, if fasta_file_object, sequences_type, infer_type or parse_method are of the wrong type.
        When calling __init__, if fasta_file_object is not a file object or is closed.
    """
    _PARSE_METHODS = ('rich', 'quick')

    def __init__(self, fasta_file_object, sequences_type=None, infer_type=False, parse_method='rich'):
        """
        Initializes file object (checks if fasta_file_object is an opened file object).

        Parameters
        ----------
        fasta_file_object : file object
            An opened file handle.
        sequences_type : 'nucleotide', 'aminoacid' or None, optional
            Indicates the type of sequences to expect ('aminoacid' or 'nucleotide'). None if unknown.
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
            If fasta_file_object, sequences_type, infer_type or parse_method are of the wrong type.
            If fasta_file_object is not a file object or is closed.
        """
        # for 'quick' parse method
        self._fasta_sequence = namedtuple('Fasta', ['header', 'sequence'])

        if hasattr(fasta_file_object, "readline"):  # assume it's a file object
            if fasta_file_object.closed and not fasta_file_object.readable():
                raise TypeError('fasta_file_object must be opened for reading')
            else:
                self._fasta_file = fasta_file_object
        else:
            raise TypeError('fasta_file_object must be a file object')

        if (isinstance(sequences_type, str) and sequences_type in LETTER_CODES) or sequences_type is None:
            self._sequences_type = sequences_type
        else:
            raise TypeError('sequence_type must be one of: %s or None' % LETTER_CODES)

        if isinstance(infer_type, bool):
            self._infer_type = infer_type
        else:
            raise TypeError('infer_type must be bool')

        if isinstance(parse_method, str) and parse_method in self._PARSE_METHODS:
            self._parse_method = parse_method
        else:
            raise TypeError('parse_method must be one of: %s' % ', '.join(self._PARSE_METHODS))

    @property
    def fasta_file(self):
        return self._fasta_file

    @property
    def sequences_type(self):
        return self._sequences_type

    @property
    def parse_method(self):
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
        if self._parse_method == 'rich':
            id_, description = self._parse_definition_line(definition_line)
            fasta_sequence = FastaSequence(sequence, id_, description, self._sequences_type, self._infer_type)
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

        definition_line = ''
        sequence = ''

        parsing_fasta_sequence = False
        for line in fasta_file:
            line = line.strip()

            # searching for '>' character at the start of a line
            # following lines contain the sequence, so parsing_fasta_sequence becomes True
            if not parsing_fasta_sequence:
                if line.startswith('>'):
                    definition_line = line
                    parsing_fasta_sequence = True

            # in the middle of parsing a FASTA sequence
            else:
                # builds sequence string until a '>' is found.
                # also ignores empty lines, which is not specified in the FASTA specification,
                # but is forgiven to badly constructed FASTA files
                if len(line) > 0 and line[0] != '>':
                    sequence += line
                elif len(line) > 0 and line[0] == '>':
                    fasta_sequence = self._generate_fasta_sequence_object(sequence, definition_line)
                    yield fasta_sequence

                    # restart variables
                    definition_line = line
                    sequence = ''
                    parsing_fasta_sequence = True

        # end of file, therefore yield last FASTA sequence
        if len(sequence) > 0:  # a FASTA sequence was actually parsed and were not just blank lines
            fasta_sequence = self._generate_fasta_sequence_object(sequence, definition_line)
            yield fasta_sequence

    def __iter__(self):
        """
        Iterates over the FASTA file.
        Returns a new iterator of the file (from the beginning) every time __iter__ is called.
        """
        if self._fasta_file.closed and not self._fasta_file.readable():  # check if file is closed
            raise TypeError('fasta_file_object must be opened for reading')

        self._current_iterator = self._iter_fasta_file(self._fasta_file)
        return self._current_iterator

    def __next__(self):
        """
        Returns the next FASTA sequence from the current iterator (most recent iterator).
        If no iterator still exists, calls __iter__ to create it.
        """
        if not hasattr(self, '_current_iterator'):
            self.__iter__()
        return next(self._current_iterator)

    def __repr__(self):
        return 'FastaParser.Reader(%s)' % self._fasta_file.name


class Writer(ParseDefinitionLine):
    """
    Writer for the given FASTA file.

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
        When calling __init__, if fasta_file_object is of the wrong type.
        When calling __init__, if fasta_file_object is not a file object or is closed.
    """

    def __init__(self, fasta_file_object):
        """
        Initializes file object (checks if fasta_file_object is a file object opened for writing).

        Parameters
        ----------
        fasta_file_object : file object
            An opened file handle ready for writing.

        Raises
        ------
        TypeError
            If fasta_file_object is of the wrong type.
            If fasta_file_object is not a file object or is closed.
        """
        if hasattr(fasta_file_object, "writelines"):  # assume it's a file object
            if fasta_file_object.closed and not fasta_file_object.writable():
                raise TypeError('fasta_file_object must be opened for writing')
            else:
                self._fasta_file = fasta_file_object
        else:
            raise TypeError('fasta_file_object must be a file object')

    @property
    def fasta_file(self):
        return self._fasta_file

    def writefasta(self, fasta_sequence):
        """
        Writes a single FASTA sequence to the provided file.

        Parameters
        ----------
        fasta_sequence : FastaSequence or (header : str, sequence : str)
            A FASTA sequence is built from the data contained in the provided FastaSequence object or the tuple of
            header + sequence.
            header may contain or not the starting '>'. header can be an empty string.
            sequence must be a non empty string.
        """
        # either use the FastaSequence object directly
        if isinstance(fasta_sequence, FastaSequence):
            pass

        # or create one with the provided header and sequence
        elif (isinstance(fasta_sequence, (tuple, list))
              and len(fasta_sequence) == 2
              and isinstance(fasta_sequence[0], str)
              and isinstance(fasta_sequence[1], str)):
            id_, description = self._parse_definition_line(fasta_sequence[0])
            sequence = ''.join(fasta_sequence[1].split('\n'))  # remove '\n's from sequence
            fasta_sequence = FastaSequence(sequence, id_, description)

        else:
            raise TypeError('fasta_sequence must be a FastaSequence object or a tuple (header : str, sequence : str)')

        # write fasta to file
        self._fasta_file.write(fasta_sequence.formatted_fasta() + '\n')

    def writefastas(self, fasta_sequences):
        """
        Writes multiple FASTA sequences to the provided file.
        Simply calls writefasta() for each object in fasta_sequences.

        Parameters
        ----------
        fasta_sequences : list of FastaSequence or list of (header : str, sequence : str)
            FASTA sequences are built from the data contained in the provided FastaSequence objects or the tuples of
            header + sequence.
            headers may contain or not the starting '>'. headers can be empty strings.
            sequences must be non empty strings.
        """
        if isinstance(fasta_sequences, (tuple, list)):
            for fasta in fasta_sequences:
                self.writefasta(fasta)
        else:
            raise TypeError('fasta_sequences must be a list of FastaSequence '
                            'objects or a list of tuples (header : str, sequence : str)')

    def __repr__(self):
        return 'FastaParser.Writer(%s)' % self._fasta_file.name
