#!python3
# coding: utf-8

"""
PyFastaParser

Parses FASTA files with the Reader class (generates FastaSequence objects or plane strings of the parsed sequences)
Writes FASTA files with the Writer class (takes FastaSequence objects or headers + sequences as strings)

ex:
    > import FastaParser
    > reader = FastaParser.Reader("fasta_file.fasta")
    > [seq.id for seq in reader]
    ['HSBGPG', 'HSGLTH1']

ex:
    > import FastaParser
    > writer = FastaParser.Writer("fasta_file.fasta")
    > seqs = [('HSBGPG example sequence', 'TTCCAGGTGTGCCAATCCAGTCCATG'),
    > ...     ('HSGLTH1 example sequence 2', 'GTACCTGACCTAACCGTGTGGACCTT')]
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
nucleotide_letter_codes_good = {
    'A': 'adenosine',
    'C': 'cytidine',
    'G': 'guanine',
    'T': 'thymidine',
    'N': 'any (A/G/C/T)',
    'U': 'uridine'
}
nucleotide_letter_codes_degenerate = {
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
nucleotide_letter_codes_complement = {
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
aminoacid_letter_codes_good = {
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
aminoacid_letter_codes_degenerate = {
    '-': 'gap of indeterminate length'
}

# set operations
letter_codes_all = set(list(nucleotide_letter_codes_good) +
                       list(nucleotide_letter_codes_degenerate) +
                       list(aminoacid_letter_codes_good) +
                       list(aminoacid_letter_codes_degenerate)
                       )
nucleotide_letter_codes_all = set(list(nucleotide_letter_codes_good) +
                                  list(nucleotide_letter_codes_degenerate))
aminoacid_letter_codes_all = set(list(aminoacid_letter_codes_good) +
                                 list(aminoacid_letter_codes_degenerate))
aminoacids_not_in_nucleotides = aminoacid_letter_codes_all - nucleotide_letter_codes_all
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
        not defined in the FASTA specification).
    supported : bool
        Indicates if letter code is supported or not
        (ie, if sequence_type is provided and letter code is defined in the FASTA specification).

    Methods
    -------
    from_lettercode()
        Alternate __init__ method. Initializes instance with a LetterCode object as only parameter.
    complement()
        Returns the complementary LetterCode of a nucleotide.

    Raises
    ------
    TypeError
        If letter_code or sequence_type are of the wrong type.
        If sequence_type is not 'nucleotide' when calling complement().
    """
    _letter_code_dictionary = {
        'nucleotide': (nucleotide_letter_codes_good, nucleotide_letter_codes_degenerate),
        'aminoacid': (aminoacid_letter_codes_good, aminoacid_letter_codes_degenerate)
    }

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
            if self._letter_code not in letter_codes_all:  # not defined in the FASTA specification
                warnings.warn('\'%s\' is not a valid letter code' % self._letter_code)
        else:
            raise TypeError('letter_code must be str of length 1')

        self._sequence_type = sequence_type
        self._description = ''
        self._degenerate = None

        if self._sequence_type in self._letter_code_dictionary:
            self._supported = True

            # letter_codes_good
            if self._letter_code in self._letter_code_dictionary[self._sequence_type][0]:
                self._description = self._letter_code_dictionary[self._sequence_type][0][self._letter_code]
                self._degenerate = False
            # letter_codes_degenerate
            elif self._letter_code in self._letter_code_dictionary[self._sequence_type][1]:
                self._description = self._letter_code_dictionary[self._sequence_type][1][self._letter_code]
                self._degenerate = True
            # _letter_code isn't defined in the FASTA specification
            else:
                self._supported = False
        elif self._sequence_type is None:
            self._supported = False
        else:
            raise TypeError('sequence_type, if defined, must be one of: %s' % ', '.join(self._letter_code_dictionary))

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

    @property
    def description(self):
        return self._description

    @property
    def degenerate(self):
        return self._degenerate

    @property
    def supported(self):
        return self._supported

    def complement(self):
        """
        Complement of nucleotide letter code.

        Returns
        -------
        LetterCode
            Complement of the current nucleotide letter code

        Raises
        ------
        TypeError
            If self.sequence_type is not 'nucleotide'.
        """
        if self._sequence_type != 'nucleotide':
            raise TypeError('Complement only works if sequence_type is \'nucleotide\'')
        return LetterCode(nucleotide_letter_codes_complement[self._letter_code], self._sequence_type)

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
        return self._letter_code


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
    complement()
        Returns the complementary FastaSequence of a nucleotide sequence.
    formatted_definition_line()
        Returns a formatted FASTA definition line (header).
    formatted_sequence(max_characters_per_line=70)
        Returns a formatted FASTA sequence (only the sequence, without the definition line).
    formatted_fasta()
        Returns a formatted FASTA (definition line and sequence).
    sequence_as_string()
        Returns the sequence as string.

    Raises
    ------
    TypeError
        If definition_line, sequence, sequence_type or infer_type are of the wrong type.
        If sequence_type is not 'nucleotide' when calling complement().
        If max_characters_per_line is not an int when calling formatted_sequence().
    """

    def __init__(self, definition_line, sequence, sequence_type=None, infer_type=False):
        """
        Initializes FASTA sequence.

        Parameters
        ----------
        definition_line: str
            FASTA sequence definition line (header). May contain or not the '>' symbol at the start.
            Can be an empty string.
        sequence : str
            String of characters representing a DNA, RNA or aminoacid sequence.
            Must be provided and cannot be empty.
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
            If definition_line, sequence, sequence_type or infer_type are of the wrong type.
        """
        if isinstance(definition_line, str):  # '>id|more_id description ...' with or without the '>' at the start
            id_and_description = definition_line.split(maxsplit=1)  # first space separates id from description

            # both id and description can be empty
            if len(id_and_description) == 0 or (len(id_and_description) == 1 and id_and_description[0] == '>'):
                self._id = ''
                self._description = ''
            else:
                if id_and_description[0].startswith('>'):
                    id_and_description[0] = id_and_description[0][1:]
                if len(id_and_description) == 1:  # description can be empty
                    self._id = id_and_description[0]
                    self._description = ''
                else:
                    self._id, self._description = id_and_description
        else:
            raise TypeError('definition_line must be str')

        if (isinstance(sequence_type, str) and sequence_type in LetterCode._letter_code_dictionary) \
                or sequence_type is None:
            self._sequence_type = sequence_type
        else:
            raise TypeError('sequence_type must be one of: %s or None' % LetterCode._letter_code_dictionary)

        if isinstance(sequence, str) and len(sequence) > 0:
            if isinstance(infer_type, bool):
                if infer_type:
                    self._sequence_type = self._infer_sequence_type(sequence)
                else:
                    self._inferred_type = False
            else:
                raise TypeError('infer_type must be bool')

            self._sequence = self._build_letter_code_sequence(sequence)
        else:
            raise TypeError('sequence must be a non empty str')

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

    @property
    def inferred_type(self):
        return self._inferred_type

    def complement(self):
        """
        Complement of nucleotide sequence.

        Returns
        -------
        FastaSequence
            Complement of the current nucleotide sequence

        Raises
        ------
        TypeError
            If self.sequence_type is not 'nucleotide'.
        """
        if self._sequence_type != 'nucleotide':
            raise TypeError('Complement only works if sequence_type is \'nucleotide\'')
        complement_sequence = ''.join([letter.complement().letter_code for letter in self])

        return FastaSequence(self._id + ' ' + self._description,
                             complement_sequence,
                             self._sequence_type)

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
            max_characters_per_line = 1 if max_characters_per_line == 0 else max_characters_per_line  # if 0

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

    def _build_letter_code_sequence(self, string_sequence):
        """
        Iterate over the given sequence and build a list of LetterCode objects.

        Parameters
        ----------
        string_sequence: str
            String of characters representing a DNA, RNA or aminoacid sequence.

        Returns
        -------
        list of LetterCode
        """
        return [LetterCode(letter_code, self._sequence_type) for letter_code in string_sequence]

    def _infer_sequence_type(self, string_sequence):
        """
        Tries to infer aminoacid sequence type.
        Tests for the presence of letter codes that can only represent aminoacids
        (ie, aminoacids letter codes not in nucleotides letter codes).
        Testing for nucleotides is not 100% accurate because there are no letter codes
        which belong solely to nucleotide type sequences.

        Parameters
        ----------
        string_sequence: str
            String of characters representing a DNA, RNA or aminoacid sequence.

        Returns
        -------
        'aminoacid' or None
            Inferred type 'aminoacid', None otherwise
        """
        for letter_code in string_sequence:
            if letter_code in aminoacids_not_in_nucleotides:
                self._inferred_type = True
                return 'aminoacid'
        self._inferred_type = False

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

    def __next__(self):
        """
        Returns the next letter code from the current iterator (most recent iterator).
        If no iterator still exists, calls __iter__ to create it.
        """
        if not hasattr(self, '_current_iterator'):
            self.__iter__()
        return next(self._current_iterator)

    def __getitem__(self, item):
        if isinstance(item, (int, slice)):
            return self._sequence[item]
        else:
            raise TypeError('Indices must be integers or slices, not tuple')

    def __len__(self):
        return len(self._sequence)

    def __repr__(self):
        """
        >id description
        sequence
        """
        return ">%s %s\n" % (self._id, self._description) + ''.join(map(str, self._sequence))

    # TODO method to return counts of each letter code in the sequence.
    #  need to implement __eq__ on LetterCode
    #  generate dict of:
    #  KEY: letter_code_string (not LetterCode objects)
    #  VALUE: counts in _build_letter_code_sequence
    #  maybe in _build_letter_code_sequence method

    # TODO method to calculate G/C content
    #  (G+C) / (A+T+G+C) * 100

    # TODO method to calculate AT/GC ratio
    #  (A+T) / (G+C)

    # TODO method to count degenerate letter codes

    # TODO Identify FASTA ID's (see linked sources)

    # TODO A new FastaSequence object must be able to be initialized with another FastaSequence object
    # TODO __get__item should return a new FastaSequence object with the sequence being the sliced result
    # TODO a slice that returns an empty sequence should return an empty FastaSequence (make sure that is possible)


class Reader:
    """
    Parser/Reader for the given FASTA file.
    Parsing mechanisms include 'quick' or 'rich':
        - 'quick': Object containing just the FASTA header and sequence attributes
        - 'rich': FastaSequence (default)

    Attributes
    ----------
    fasta_file : FASTA file object
        The FASTA file passed as parameter.
    sequences_type : 'nucleotide', 'aminoacid' or None
        Indicates the type of sequences to expect ('aminoacid' or 'nucleotide'). Can be None if not known.
    infer_type: bool
        True if sequence type of each sequence is to be inferred, False otherwise.
    parse_method: 'rich' or 'quick'
        Parse method used ('rich' or 'quick').

    Raises
    ------
    TypeError
        If fasta_file_object, sequences_type, infer_type or parse_method are of the wrong type.
        If fasta_file_object is not a file object or is closed.
    """

    def __init__(self, fasta_file_object, sequences_type=None, infer_type=False, parse_method='rich'):
        """
        Initializes file object (checks if fasta_file_object is an opened file object).
        Only one FASTA sequence is read at a time from the FASTA file, previous sequences are not kept in memory.

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
            infer_type parameters are ignored
            'rich' implements more functionality, but is a bit slower.

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

        if (isinstance(sequences_type, str) and sequences_type in LetterCode._letter_code_dictionary) \
                or sequences_type is None:
            self._sequences_type = sequences_type
        else:
            raise TypeError('sequence_type must be one of: %s or None' % LetterCode._letter_code_dictionary)

        if isinstance(infer_type, bool):
            self._infer_type = infer_type
        else:
            raise TypeError('infer_type must be bool')

        parse_methods = ('rich', 'quick')
        if isinstance(parse_method, str) and parse_method in parse_methods:
            self._parse_method = parse_method
        else:
            raise TypeError('parse_method must be one of: %s' % ', '.join(parse_methods))

    @property
    def fasta_file(self):
        return self._fasta_file

    @property
    def sequences_type(self):
        return self._sequences_type

    @property
    def infer_type(self):
        return self._infer_type

    @property
    def parse_method(self):
        return self._parse_method

    def __iter__(self):
        """
        Iterates over the FASTA file.
        Returns a new iterator of the file (from the beginning) every time __iter__ is called.
        """
        if self._fasta_file.closed:  # check if file is closed
            raise TypeError('fasta_file_object must be opened for reading')

        def iter_fasta_file(fasta_file):
            fasta_file.seek(0)  # restart cursor position (just in case)

            definition_line = ''
            sequence = ''

            parsing_fasta_sequence = False
            for line in fasta_file:
                line = line.strip()

                # searching for '>' character at the start of a line
                if not parsing_fasta_sequence:
                    if line.startswith('>'):
                        definition_line = line
                        parsing_fasta_sequence = True

                # in the middle of parsing a FASTA sequence
                else:
                    if len(line) > 0 and line[0] != '>':
                        sequence += line
                    elif len(line) > 0 and line[0] == '>':
                        if self._parse_method == 'rich':
                            fasta_sequence = FastaSequence(definition_line, sequence, self._sequences_type,
                                                           self._infer_type)
                        else:  # 'quick'
                            fasta_sequence = self._fasta_sequence(definition_line, sequence)
                        yield fasta_sequence

                        # restart variables
                        definition_line = line
                        sequence = ''
                        parsing_fasta_sequence = True

            # end of file, therefore yield last FASTA sequence
            if len(sequence) > 0:  # a FASTA sequence was actually parsed and were not just blank lines
                if self._parse_method == 'rich':
                    fasta_sequence = FastaSequence(definition_line, sequence, self._sequences_type, self._infer_type)
                else:  # 'quick'
                    fasta_sequence = self._fasta_sequence(definition_line, sequence)
                yield fasta_sequence

        self._current_iterator = iter_fasta_file(self._fasta_file)
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
        return "<%s - FASTAFILE:%s>" % (self.__class__.__name__, self._fasta_file.name)


class Writer:
    """
    Writer for the given FASTA file.

    Attributes
    ----------
    fasta_file : FASTA file object
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
        If fasta_file_object is of the wrong type.
        If fasta_file_object is not a file object or is closed.
    """

    def __init__(self, fasta_file_object):
        """
        Initializes file object (checks if fasta_file_object is a file object opened for writing).

        Parameters
        ----------
        fasta_file_object : file object
            An opened file handle.

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
        elif isinstance(fasta_sequence, (tuple, list)) and len(fasta_sequence) == 2 \
                and isinstance(fasta_sequence[0], str) and isinstance(fasta_sequence[1], str):
            header = fasta_sequence[0]
            sequence = ''.join(fasta_sequence[1].split('\n'))  # remove '\n's from sequence
            fasta_sequence = FastaSequence(header, sequence)

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
