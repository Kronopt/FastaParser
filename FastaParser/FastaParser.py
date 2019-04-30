#!python3
# coding: utf-8

"""
PyFastaParser

Parses FASTA files with the FastaParser class and generates
FastaSequence objects of the parsed sequences

ex:
    > import FastaParser
    > parser = FastaParser.FastaParser("fasta_file.fasta")
    > [seq.id for seq in parser]
    ['HSBGPG', 'HSGLTH1']


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
        if isinstance(letter_code, str):
            self._letter_code = letter_code.upper()
            if self._letter_code not in letter_codes_all:  # not defined in the FASTA specification
                warnings.warn('\'%s\' is not a valid letter code' % self._letter_code)
        else:
            raise TypeError('letter_code must be str')
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

    def __repr__(self):
        return self._letter_code


class FastaSequence:
    """
    Represents one FASTA sequence.
    The class itself is an iterator of the sequence it represents.

    Attributes
    ----------
    id : str
        ID portion of the definition line (header).
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

    Raises
    ------
    TypeError
        If definition_line, sequence, sequence_type or infer_type are of the wrong type.
        If sequence_type is not 'nucleotide' when calling complement().
    """

    def __init__(self, definition_line, sequence, sequence_type=None, infer_type=False):
        """
        Initializes FASTA sequence.

        Parameters
        ----------
        definition_line: str
            FASTA sequence definition line (header), containing the '>' symbol at the start.
        sequence : str
            String of characters representing a DNA, RNA or aminoacid sequence.
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
        if isinstance(definition_line, str):  # '>id|more_id description ...'
            id_and_description = definition_line.split(maxsplit=1)  # first space separates id from description
            if len(id_and_description) == 1 and id_and_description[0] == '>':  # both id and description can be empty
                self._id = ''
                self._description = ''
            else:
                # assumes a starting '>'
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

        if isinstance(sequence, str):
            if isinstance(infer_type, bool):
                if infer_type:
                    self._sequence_type = self._infer_sequence_type(sequence)
                else:
                    self._inferred_type = False
            else:
                raise TypeError('infer_type must be bool')

            self._sequence = self._build_letter_code_sequence(sequence)
        else:
            raise TypeError('sequence must be a str')

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

        return FastaSequence('>' + self._id + ' ' + self._description,
                             complement_sequence,
                             self._sequence_type)

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
        """
        def iter_sequence():
            for letter_code in self._sequence:
                yield letter_code
        return iter_sequence()

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

    # TODO method to return a string containing a formatted FASTA sequence
    #  with header and sequence lines (maximum of 80 characters per sequence line)

    # TODO document non protected methods (if any)


class FastaParser:
    """
    Parses the given FASTA file.
    Parsing methods include 'quick' or 'rich':
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
            Indicates if FastaParser should try to infer aminoacid sequence type for each sequence.
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
            if fasta_file_object.closed:
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

        return iter_fasta_file(self._fasta_file)

    def __repr__(self):
        return "<%s - FASTAFILE:%s>" % (self.__class__.__name__, self._fasta_file.name)

    # TODO FASTA writer
    # TODO Identify ID's (see linked sources)
