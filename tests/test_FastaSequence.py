#!python
# coding: utf-8

"""
Tests for FastaParser.FastaSequence class
"""


import pytest
from FastaParser import FastaSequence, LetterCode, \
    NUCLEOTIDE_LETTER_CODES_GOOD, AMINOACID_LETTER_CODES_GOOD, \
    NUCLEOTIDE_LETTER_CODES_DEGENERATE, AMINOACID_LETTER_CODES_DEGENERATE, \
    NUCLEOTIDE_LETTER_CODES_COMPLEMENT, AMINOACIDS_NOT_IN_NUCLEOTIDES


##########
# Fixtures
##########


@pytest.fixture()
def nucleotide_good():
    return (FastaSequence(''.join(NUCLEOTIDE_LETTER_CODES_GOOD), sequence_type='nucleotide'),
            [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD])


@pytest.fixture()
def nucleotide_degenerate():
    return (FastaSequence(''.join(NUCLEOTIDE_LETTER_CODES_DEGENERATE), sequence_type='nucleotide'),
            [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_DEGENERATE])


@pytest.fixture()
def nucleotide_good_complement():
    return [LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT[letter_code]) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD]


@pytest.fixture()
def aminoacid_good():
    return (FastaSequence(''.join(AMINOACID_LETTER_CODES_GOOD), sequence_type='aminoacid'),
            [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_GOOD])


@pytest.fixture()
def aminoacid_degenerate():
    return (FastaSequence(''.join(AMINOACID_LETTER_CODES_DEGENERATE), sequence_type='aminoacid'),
            [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_DEGENERATE])


@pytest.fixture()
def sequence_type_none():
    return (FastaSequence(''.join(NUCLEOTIDE_LETTER_CODES_GOOD) + ''.join(AMINOACID_LETTER_CODES_GOOD)),
            [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD] +
            [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_GOOD])


@pytest.fixture()
def sequence_type_none_complement():
    return ([LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT[letter_code])
            for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD] +
            [LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT.get(letter_code, letter_code))
             for letter_code in AMINOACID_LETTER_CODES_GOOD])


@pytest.fixture()
def letter_codes_unknown(unknown_characters):
    return (FastaSequence(''.join(unknown_characters)),
            [LetterCode(letter_code) for letter_code in unknown_characters])


@pytest.fixture()
def actg_letter_code_list():
    return [LetterCode(letter_code) for letter_code in 'ACTG']


#######
# Tests
#######


class Test__Init__:
    def test_sequence_nucleotide_good(self, nucleotide_good):
        fasta_sequence, correct_sequence = nucleotide_good
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'nucleotide'
        assert fasta_sequence.inferred_type is False

    def test_sequence_nucleotide_degenerate(self, nucleotide_degenerate):
        fasta_sequence, correct_sequence = nucleotide_degenerate
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'nucleotide'
        assert fasta_sequence.inferred_type is False

    def test_sequence_aminoacid_good(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'aminoacid'
        assert fasta_sequence.inferred_type is False

    def test_sequence_aminoacid_degenerate(self, aminoacid_degenerate):
        fasta_sequence, correct_sequence = aminoacid_degenerate
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'aminoacid'
        assert fasta_sequence.inferred_type is False

    def test_sequence_unknown_characters(self, letter_codes_unknown):
        fasta_sequence, correct_sequence = letter_codes_unknown
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_sequence_not_str(self):
        with pytest.raises(TypeError):
            FastaSequence(1)
        with pytest.raises(TypeError):
            FastaSequence([])
        with pytest.raises(TypeError):
            FastaSequence(FastaSequence('actg'))

    def test_sequence_empty_str(self):
        with pytest.raises(TypeError):
            FastaSequence('')

    def test_sequence_lower_case(self):
        fasta_sequence = FastaSequence(''.join(NUCLEOTIDE_LETTER_CODES_GOOD).lower())
        correct_sequence = [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_id_good(self, actg_letter_code_list):
        id_ = 'correct_id|some_other_id'
        # include '>'
        fasta_sequence = FastaSequence('ACTG', id_='>' + id_)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # exclude '>'
        fasta_sequence = FastaSequence('ACTG', id_=id_)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_id_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence('ACTG', id_=1)
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence('ACTG', id_=[])

    # def test_id_empty_str (already tested)

    def test_description_good(self, actg_letter_code_list):
        description = 'some correct description with numbers 1234567890#'
        fasta_sequence = FastaSequence('ACTG', description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_description_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence('ACTG', description=1)
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence('ACTG', description=[])

    # def test_description_empty_str (already tested)

    def test_id_and_description_good(self, actg_letter_code_list):
        id_ = 'correct_id|some_other_id'
        description = 'some correct description with numbers 1234567890#'
        # include '>'
        fasta_sequence = FastaSequence('ACTG', id_='>' + id_, description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # exclude '>'
        fasta_sequence = FastaSequence('ACTG', id_=id_, description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_sequence_type_nucleotide(self, unknown_characters, letter_codes_unknown):
        # nucleotide_good already tested
        fasta_sequence = FastaSequence(''.join(unknown_characters), sequence_type='nucleotide')
        assert fasta_sequence.sequence == letter_codes_unknown[1]
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'nucleotide'
        assert fasta_sequence.inferred_type is False

    def test_sequence_type_aminoacid(self, unknown_characters, letter_codes_unknown):
        # aminoacid_good already tested
        fasta_sequence = FastaSequence(''.join(unknown_characters), sequence_type='aminoacid')
        assert fasta_sequence.sequence == letter_codes_unknown[1]
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'aminoacid'
        assert fasta_sequence.inferred_type is False

    # def test_sequence_type_none (already tested)

    def test_sequence_type_wrong_type(self):
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type=1)
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type=[])

    def test_sequence_type_wrong_str(self):
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type='')
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type='wrong string')
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type='aminoacids')
        with pytest.raises(TypeError):
            FastaSequence('ACTG', sequence_type='nucleotides')

    def test_infer_type_true(self, nucleotide_good, aminoacid_good):
        # nucleotide sequence
        fasta_sequence = FastaSequence(''.join(NUCLEOTIDE_LETTER_CODES_GOOD), infer_type=True)
        correct_sequence = nucleotide_good[1]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # aminoacid sequence
        fasta_sequence = FastaSequence(''.join(AMINOACID_LETTER_CODES_GOOD), infer_type=True)
        correct_sequence = aminoacid_good[1]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'aminoacid'
        assert fasta_sequence.inferred_type is True

    # def test_infer_type_false (already tested)

    def test_infer_type_not_bool(self):
        with pytest.raises(TypeError):
            FastaSequence('ACTG', infer_type=1)
        with pytest.raises(TypeError):
            FastaSequence('ACTG', infer_type=[])
        with pytest.raises(TypeError):
            FastaSequence('ACTG', infer_type='')

    def test_counts_good(self, nucleotide_good, aminoacid_good):
        assert nucleotide_good[0]._counts == dict(zip(NUCLEOTIDE_LETTER_CODES_GOOD,
                                                      [1]*len(NUCLEOTIDE_LETTER_CODES_GOOD)))
        assert aminoacid_good[0]._counts == dict(zip(AMINOACID_LETTER_CODES_GOOD,
                                                     [1]*len(AMINOACID_LETTER_CODES_GOOD)))

    def test_counts_same_character(self):
        sequence = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        fasta_sequence = FastaSequence(sequence)
        assert fasta_sequence._counts == {'A': len(sequence)}

    def test_counts_letter_codes_unknown(self, letter_codes_unknown, unknown_characters):
        assert letter_codes_unknown[0]._counts == dict(zip(unknown_characters, [1]*len(unknown_characters)))

    def test_gc_content(self, nucleotide_good):
        assert nucleotide_good[0]._gc_content is None

    def test_at_gc_ratio(self, nucleotide_good):
        assert nucleotide_good[0]._at_gc_ratio is None


class Test_from_fastasequence:
    def test_fastasequence_good(self, nucleotide_good):
        original = nucleotide_good[0]
        new = FastaSequence.from_fastasequence(original)
        assert new is not original
        assert new.sequence == original.sequence
        assert new.id == original.id
        assert new.description == original.description
        assert new.sequence_type == original.sequence_type
        assert new.inferred_type == original.inferred_type

    def test_fastasequence_wrong_type(self):
        with pytest.raises(TypeError):
            FastaSequence.from_fastasequence('ACTG')
        with pytest.raises(TypeError):
            FastaSequence.from_fastasequence(1)


class Test_sequence_type_property:
    def test_set_nucleotide(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        fasta_sequence.sequence_type = 'nucleotide'
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'nucleotide'
        assert fasta_sequence.inferred_type is False

    # def test_get_nucleotide (already tested in Test__Init__)

    def test_set_aminoacid(self, nucleotide_good):
        fasta_sequence, correct_sequence = nucleotide_good
        fasta_sequence.sequence_type = 'aminoacid'
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type == 'aminoacid'
        assert fasta_sequence.inferred_type is False

    # def test_get_aminoacid (already tested in Test__Init__)

    def test_set_none(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        fasta_sequence.sequence_type = None
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ''
        assert fasta_sequence.description == ''
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_none (already tested in Test__Init__)

    def test_set_wrong_str(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = 'something'
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = ''

    # def test_get_wrong_str (already tested in Test__Init__)

    def test_set_wrong_type(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = 1
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = []

    # def test_get_wrong_type (already tested in Test__Init__)


class Test_complement:
    def test_sequence_type_aminoacid(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.complement()

    def test_sequence_type_nucleotide_reverse_False(self, nucleotide_good, nucleotide_good_complement):
        fasta_sequence = nucleotide_good[0]
        complement = fasta_sequence.complement()
        assert complement.sequence == nucleotide_good_complement
        assert complement.id == fasta_sequence.id
        assert complement.description == '[COMPLEMENT]'
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_sequence_type_nucleotide_reverse_true(self, nucleotide_good, nucleotide_good_complement):
        fasta_sequence = nucleotide_good[0]
        complement = fasta_sequence.complement(True)
        assert complement.sequence == nucleotide_good_complement[::-1]
        assert complement.id == fasta_sequence.id
        assert complement.description == '[REVERSE COMPLEMENT]'
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_reverse_not_bool(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.complement(1)
        with pytest.raises(TypeError):
            fasta_sequence.complement([])
        with pytest.raises(TypeError):
            fasta_sequence.complement('')

    def test_sequence_type_none(self, sequence_type_none, sequence_type_none_complement):
        fasta_sequence = sequence_type_none[0]
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == sequence_type_none_complement
        assert complement.id == fasta_sequence.id
        assert complement.description == '[COMPLEMENT]'
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type
        # aminoacid letter codes that are not also a nucleotide
        fasta_sequence = FastaSequence(''.join(AMINOACIDS_NOT_IN_NUCLEOTIDES))
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == [LetterCode(letter_code) for letter_code in AMINOACIDS_NOT_IN_NUCLEOTIDES]
        assert complement.id == fasta_sequence.id
        assert complement.description == '[COMPLEMENT]'
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_letter_code_unknown(self, letter_codes_unknown):
        fasta_sequence, correct_sequence = letter_codes_unknown
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == correct_sequence
        assert complement.id == fasta_sequence.id
        assert complement.description == '[COMPLEMENT]'
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type


# tested in Test__Init__:
#   class Test_id_property
#   class Test_description_property
#   class Test_sequence_property
#   class Test_inferred_type


# TODO
# class Test_gc_content
# class Test_at_gc_ratio
# class Test_count_letter_codes
# class Test_count_letter_codes_degenerate
# class Test_formatted_definition_line
# class Test_formatted_sequence
# class Test_formatted_fasta
# class Test_sequence_as_string
# class Test_reverse
# class Test__iter__
# class Test__reversed__
# class Test__next__
# class Test__getitem__
# class Test__len__
# class Test__repr__
# class Test__str__
