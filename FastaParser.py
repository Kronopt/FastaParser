#!python3
# coding: utf-8

"""
PyFastaParser

Parses fasta files with FastaParser class and generates FastaSequence objects of the parsed sequences

ex:
    > import FastaParser
    > parser = FastaParser.FastaParser("fasta_file.fasta")
    > [seq.id for seq in parser]
    ['HSBGPG', 'HSGLTH1']
"""

__author__ = 'Pedro HC David, https://github.com/Kronopt'
__credits__ = ['Pedro HC David']
__version__ = '1.0'


class FastaSequence:
    """
    Represents one sequence as read from the given fasta file.
    The class itself is an iterator of the sequence it contains.
    """

    def __init__(self, seq_id, description, sequence):
        """
        Simply initializes the id, description and sequence.

        Parameters
        ----------
        seq_id : str
            ID portion of the sequence.
        description : str
            Description portion of the string. Can be empty.
        sequence : str
            Sequence.
        """
        self.id = seq_id
        self.description = description
        self.sequence = sequence

    def __iter__(self):
        """
        Iterates over the sequence.
        """
        def iter_sequence(sequence):
            for character in sequence:
                yield character
        return iter_sequence(self.sequence)

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "<%s object - ID:%s | DESCRIPTION:%s>" % (self.__class__.__name__, self.id, self.description)


class FastaParser:
    """
    Parses the given fasta file.
    Assumes fasta file is properly formatted.
    """

    def __init__(self, fasta_file, keep_sequences=False):
        """
        Initializes file object (checks if fasta_file is a path or a file object).
        Only one sequence at a time is read from the fasta file, previous sequences are not kept in memory.
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
        self.sequences = []

    def __iter__(self):
        """
        Iterates over the fasta file.
        """
        # check if file was closed, and open it again for new iteration
        if self._fasta.closed:
            self._fasta = open(self._fasta.name, 'r')

        def iter_fasta_file(fasta_file):
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
                        self._iteration_ended = True  # Iterated once. No more sequences will be added to self.sequences
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
                    self.sequences.append(fasta_sequence)

                yield fasta_sequence

        return iter_fasta_file(self._fasta)

    def __repr__(self):
        return "<%s object - FASTAFILE:%s>" % (self.__class__.__name__, self._fasta.name)
