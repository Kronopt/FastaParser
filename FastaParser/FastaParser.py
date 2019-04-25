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
"""

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

    Raises
    ------
    TypeError
        If letter_code is not str.
    ValueError
        If sequence_type is not 'nucleotide', 'aminoacid' or None.
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
        sequence_type : str or None, optional
            'nucleotide' or 'aminoacid.
        """
        if isinstance(letter_code, str):
            self._letter_code = letter_code.upper()
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
            raise ValueError('sequence_type, if defined, must be one of: %s' % ', '.join(self._letter_code_dictionary))

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
    sequence_type : str or None
        Indicates the type of sequence ('aminoacid' or 'nucleotide'). Can be None if not known.
    inferred_type: bool
        True if FastaSequence inferred the sequence type, False otherwise.

    Raises
    ------
    TypeError
        If definition_line, sequence, sequence_type or infer_type are of the wrong type.
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
        sequence_type : str or None, optional
            Indicates the type of sequence ('aminoacid' or 'nucleotide').
            If not defined, FastaSequence can try to infer type based on the letter codes.
        infer_type : bool, optional
            Indicates if FastaSequence should try to infer aminoacid sequence type.
            If True, FastaSequence will analyse the whole sequence, in the worst case scenario,
            and can only identify aminoacid sequences.
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
            raise TypeError('sequence_type must be str or None')

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
        str or None
            Inferred type
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
        <FastaSequence - ID:id#123 | DESCRIPTION:sequence description trimmed to 50 characters>
        """
        _id = self._id if self._id else '\'\''
        _description = self._description[:50] if self.description else '\'\''
        return "<%s - ID:%s | DESCRIPTION:%s>" % (self.__class__.__name__, _id, _description)

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

    # TODO document non protected methods


# TODO define the following in FastaParser:
# TODO sequence_type
# TODO infer_type

# TODO redo
class FastaParser:
    """
    Parses the given FASTA file.
    Assumes FASTA file is properly formatted.
    """

    def __init__(self, fasta_file, keep_sequences=False):
        """
        Initializes file object (checks if fasta_file is a path or a file object).
        Only one sequence at a time is read from the FASTA file, previous sequences are not kept in memory.
        If keep_sequences=True, previous sequences are kept in memory.

        Parameters
        ----------
        fasta_file : str / file
            path string or a file object
        keep_sequences : bool
            Whether to keep the sequences iterated over or not
        """
        if isinstance(fasta_file, str):  # try to open file path
            self._fasta = open(fasta_file, 'r')
        elif hasattr(fasta_file, "readline"):  # assume it's a file object
            if fasta_file.closed:
                self._fasta = open(fasta_file.name, 'r')
            else:
                self._fasta = fasta_file
        else:
            raise Exception("Not a file.")

        self._iteration_ended = False
        self._keep_sequences = keep_sequences
        self.sequences = {}

    def __iter__(self):
        """
        Iterates over the FASTA file.
        """
        # check if file was closed, and open it again for new iteration
        if self._fasta.closed:
            self._fasta = open(self._fasta.name, 'r')

        def iter_fasta_file(fasta_file):
            # TODO use for loop instead
            # TODO check if sequence name is unique (ie, check if sequence name is already in dict)
            # assumes first line begins with '>'
            first_line = fasta_file.readline()[1:].split(" ", 1)

            end_of_file = False
            while not end_of_file:
                if len(first_line) == 1:  # description can be empty
                    seq_id = first_line[0].rstrip()
                    seq_description = ''
                else:
                    seq_id, seq_description = first_line

                seq = ''
                end_of_sequence = False
                while not end_of_sequence:
                    sequence_line = fasta_file.readline()
                    if len(sequence_line) == 0:  # end of file, end iteration
                        end_of_sequence = True
                        end_of_file = True
                        fasta_file.close()
                    elif not sequence_line.startswith('>'):  # Line containing part of the sequence
                        seq += sequence_line.rstrip()
                    else:
                        end_of_sequence = True
                        first_line = sequence_line[1:].split(" ", 1)

                fasta_sequence = FastaSequence(seq_id, seq_description.rstrip(), seq)

                if self._keep_sequences and not self._iteration_ended:
                    self.sequences[fasta_sequence.id] = fasta_sequence

                yield fasta_sequence
            
            self._iteration_ended = True  # Iterated once. No more sequences will be added to self.sequences

        return iter_fasta_file(self._fasta)

    def __repr__(self):
        return "<%s - FASTAFILE:%s>" % (self.__class__.__name__, self._fasta.name)
