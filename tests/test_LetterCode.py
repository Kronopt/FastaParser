#!python
# coding: utf-8

"""
Tests for FastaParser.LetterCode class
"""


import pytest
from FastaParser import LetterCode


##########
# Fixtures
##########


@pytest.fixture(scope='module')
def nucleotide_good():
    return LetterCode('A', 'nucleotide')


@pytest.fixture(scope='module')
def nucleotide_degenerate():
    return LetterCode('K', 'nucleotide')


@pytest.fixture(scope='module')
def aminoacid_good():
    return LetterCode('H', 'aminoacid')


@pytest.fixture(scope='module')
def aminoacid_degenerate():
    return LetterCode('-', 'aminoacid')


@pytest.fixture(scope='module')
def sequence_type_none():
    return LetterCode('A')


@pytest.fixture(scope='module')
def unknown_characters():
    return 'O', '»', '%', 'º', '?', 'غ', '\n'


@pytest.fixture(scope='module')
def letter_codes_unknown(unknown_characters):
    letter_codes = [LetterCode(character) for character in unknown_characters]
    return zip(letter_codes, unknown_characters)


@pytest.fixture(scope='module')
def letter_codes_unknown_nucleotide(unknown_characters):
    letter_codes = [LetterCode(character, 'nucleotide') for character in unknown_characters]
    return zip(letter_codes, unknown_characters)


@pytest.fixture(scope='module')
def letter_codes_unknown_aminoacid(unknown_characters):
    letter_codes = [LetterCode(character, 'aminoacid') for character in unknown_characters]
    return zip(letter_codes, unknown_characters)


#######
# Tests
#######


class Test__Init__:
    def test_letter_code_good(self, nucleotide_good, aminoacid_good):
        # nucleotide
        assert nucleotide_good.letter_code == 'A'
        assert nucleotide_good.sequence_type == 'nucleotide'
        assert nucleotide_good.description == 'adenosine'
        assert nucleotide_good.degenerate is False
        assert nucleotide_good.supported is True
        # aminoacid
        assert aminoacid_good.letter_code == 'H'
        assert aminoacid_good.sequence_type == 'aminoacid'
        assert aminoacid_good.description == 'histidine'
        assert aminoacid_good.degenerate is False
        assert aminoacid_good.supported is True

    def test_letter_code_degenerate(self, nucleotide_degenerate, aminoacid_degenerate):
        # nucleotide
        assert nucleotide_degenerate.letter_code == 'K'
        assert nucleotide_degenerate.sequence_type == 'nucleotide'
        assert nucleotide_degenerate.description == 'keto (G/T)'
        assert nucleotide_degenerate.degenerate is True
        assert nucleotide_degenerate.supported is True
        # aminoacid
        assert aminoacid_degenerate.letter_code == '-'
        assert aminoacid_degenerate.sequence_type == 'aminoacid'
        assert aminoacid_degenerate.description == 'gap of indeterminate length'
        assert aminoacid_degenerate.degenerate is True
        assert aminoacid_degenerate.supported is True

    def test_letter_code_unknown_characters(self, letter_codes_unknown):
        for letter_code, character in letter_codes_unknown:
            assert letter_code.letter_code == character
            assert letter_code.sequence_type is None
            assert letter_code.description == ''
            assert letter_code.degenerate is None
            assert letter_code.supported is False

    def test_letter_code_not_str(self):
        with pytest.raises(TypeError):
            LetterCode(1)
        with pytest.raises(TypeError):
            LetterCode(LetterCode('A'))

    def test_letter_code_multiple_characters(self):
        with pytest.raises(TypeError):
            LetterCode('AC')

    def test_letter_code_empty_str(self):
        with pytest.raises(TypeError):
            LetterCode('')

    def test_letter_code_lower_case(self):
        letter_code = LetterCode('a', 'nucleotide')
        assert letter_code.letter_code == 'A'
        assert letter_code.sequence_type == 'nucleotide'
        assert letter_code.description == 'adenosine'
        assert letter_code.degenerate is False
        assert letter_code.supported is True

    def test_sequence_type_nucleotide(self, letter_codes_unknown_nucleotide):
        # good and degenerate already tested in test_letter_code_good and test_letter_code_degenerate
        for letter_code, character in letter_codes_unknown_nucleotide:
            assert letter_code.letter_code == character
            assert letter_code.sequence_type == 'nucleotide'
            assert letter_code.description == ''
            assert letter_code.degenerate is None
            assert letter_code.supported is False

    def test_sequence_type_aminoacid(self, letter_codes_unknown_aminoacid):
        # good and degenerate already tested in test_letter_code_good and test_letter_code_degenerate
        for letter_code, character in letter_codes_unknown_aminoacid:
            assert letter_code.letter_code == character
            assert letter_code.sequence_type == 'aminoacid'
            assert letter_code.description == ''
            assert letter_code.degenerate is None
            assert letter_code.supported is False

    def test_sequence_type_none(self, sequence_type_none):
        assert sequence_type_none.letter_code == 'A'
        assert sequence_type_none.sequence_type is None
        assert sequence_type_none.description == ''
        assert sequence_type_none.degenerate is None
        assert sequence_type_none.supported is False

    def test_sequence_type_incorrect(self):
        with pytest.raises(TypeError):
            LetterCode('A', 1)

        with pytest.raises(TypeError):
            LetterCode('A', [])

        with pytest.raises(TypeError):
            LetterCode('A', '')

        with pytest.raises(TypeError):
            LetterCode('A', 'aminoacidnucleotide')


class Test_from_lettercode:
    def test_lettercode_good(self):
        pass

    def test_lettercode_wrong_type(self):
        pass


class Test_letter_code_property:
    def test_nucleotide_good(self):
        pass

    def test_aminoacid_good(self):
        pass

    def test_unknown_characters(self):
        pass

    def test_weird_characters(self):
        # \n # $ % ? . -, unicode characters, etc etc
        pass


class Test_sequence_type_property:
    def test_get_nucleotide(self):
        pass

    def test_get_aminoacid(self):
        pass

    def test_get_none(self):
        pass

    def test_set_nucleotide(self):
        pass

    def test_set_aminoacid(self):
        pass

    def test_set_none(self):
        pass

    def test_set_wrong_str(self):
        pass

    def test_set_wrong_type(self):
        pass


class Test_description_property:
    def test_sequence_type_nucleotide_good(self):
        pass

    def test_sequence_type_nucleotide_degenerate(self):
        pass

    def test_sequence_type_aminoacid_good(self):
        pass

    def test_sequence_type_aminoacid_degenerate(self):
        pass

    def test_sequence_type_none(self):
        pass

    def test_letter_code_unknown(self):
        pass


class Test_degenerate_property:
    def test_letter_code_degenerate(self):
        pass

    def test_letter_code_good(self):
        pass

    def test_letter_code_unknown(self):
        pass

    def test_sequence_type_none(self):
        pass


class Test_supported_property:
    def test_letter_code_degenerate(self):
        pass

    def test_letter_code_good(self):
        pass

    def test_letter_code_unknown(self):
        pass

    def test_sequence_type_none(self):
        pass


class Test_complement:
    def test_sequence_type_aminoacid(self):
        # throws error
        pass

    def test_sequence_type_nucleotide(self):
        pass

    def test_sequence_type_none(self):
        pass

    def test_letter_code_unknown(self):
        # complemente == letter_code
        pass


class Test__eq__:
    def test_lettercode(self):
        pass

    def test_str(self):
        pass

    def test_not_str_or_lettercode(self):
        pass


class Test__repr__:
    def test__repr__(self):
        pass


class Test__str__:
    def test__str__(self):
        pass
