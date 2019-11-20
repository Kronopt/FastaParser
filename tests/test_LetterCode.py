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


@pytest.fixture(scope='function')
def nucleotide_good():
    return LetterCode('A', 'nucleotide')


@pytest.fixture(scope='module')
def nucleotide_degenerate():
    return LetterCode('K', 'nucleotide')


@pytest.fixture(scope='function')
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
    def test_lettercode_good(self, nucleotide_good):
        new = LetterCode.from_lettercode(nucleotide_good)
        assert new is not nucleotide_good
        assert new.letter_code == nucleotide_good.letter_code
        assert new.sequence_type == nucleotide_good.sequence_type
        assert new.description == nucleotide_good.description
        assert new.degenerate == nucleotide_good.degenerate
        assert new.supported == nucleotide_good.supported

    def test_lettercode_wrong_type(self):
        with pytest.raises(TypeError):
            LetterCode.from_lettercode('A')
        with pytest.raises(TypeError):
            LetterCode.from_lettercode(1)


class Test_sequence_type_property:
    # def test_get_nucleotide  # tested in Test__Init__
    # def test_get_aminoacid
    # def test_get_none

    def test_set_nucleotide(self, aminoacid_good):
        aminoacid_good.sequence_type = 'nucleotide'
        assert aminoacid_good.letter_code == 'H'
        assert aminoacid_good.sequence_type == 'nucleotide'
        assert aminoacid_good.description == 'A/C/T'
        assert aminoacid_good.degenerate is True
        assert aminoacid_good.supported is True

    def test_set_aminoacid(self, nucleotide_good):
        nucleotide_good.sequence_type = 'aminoacid'
        assert nucleotide_good.letter_code == 'A'
        assert nucleotide_good.sequence_type == 'aminoacid'
        assert nucleotide_good.description == 'alanine'
        assert nucleotide_good.degenerate is False
        assert nucleotide_good.supported is True

    def test_set_none(self, aminoacid_good):
        aminoacid_good.sequence_type = None
        assert aminoacid_good.letter_code == 'H'
        assert aminoacid_good.sequence_type is None
        assert aminoacid_good.description == ''
        assert aminoacid_good.degenerate is None
        assert aminoacid_good.supported is False

    def test_set_wrong_str(self, aminoacid_good):
        with pytest.raises(TypeError):
            aminoacid_good.sequence_type = 'something'
        with pytest.raises(TypeError):
            aminoacid_good.sequence_type = ''

    def test_set_wrong_type(self, aminoacid_good):
        with pytest.raises(TypeError):
            aminoacid_good.sequence_type = 1
        with pytest.raises(TypeError):
            aminoacid_good.sequence_type = []


class Test_complement:
    def test_sequence_type_aminoacid(self, aminoacid_good):
        with pytest.raises(TypeError):
            aminoacid_good.complement()

    def test_sequence_type_nucleotide(self, nucleotide_good):
        complement = nucleotide_good.complement()
        assert complement.letter_code == 'T'
        assert complement.sequence_type == 'nucleotide'
        assert complement.description == 'thymidine'
        assert complement.degenerate is False
        assert complement.supported is True

    def test_sequence_type_none(self, sequence_type_none):
        complement = sequence_type_none.complement()
        assert complement.letter_code == 'T'
        assert complement.sequence_type is None
        assert complement.description == ''
        assert complement.degenerate is None
        assert complement.supported is False
        # aminoacid letter code that is not also a nucleotide
        aminoacid = LetterCode('X').complement()
        assert aminoacid.letter_code == 'X'
        assert aminoacid.sequence_type is None
        assert aminoacid.description == ''
        assert aminoacid.degenerate is None
        assert aminoacid.supported is False

    def test_letter_code_unknown(self, letter_codes_unknown):
        for letter_code, character in letter_codes_unknown:
            complement = letter_code.complement()
            assert complement.letter_code == character
            assert complement.sequence_type is None
            assert complement.description == ''
            assert complement.degenerate is None
            assert complement.supported is False


class Test__eq__:
    def test_lettercode(self, nucleotide_good):
        assert nucleotide_good == LetterCode('A')
        assert nucleotide_good != LetterCode('C')

    def test_str(self, nucleotide_good):
        assert nucleotide_good == 'A'
        assert nucleotide_good == 'a'
        assert nucleotide_good != 'C'
        assert nucleotide_good != ''

    def test_not_str_or_lettercode(self):
        assert nucleotide_good != 1
        assert nucleotide_good != []


class Test__repr__:
    def test__repr__(self, nucleotide_good):
        assert repr(nucleotide_good) == 'LetterCode(\'A\')'


class Test__str__:
    def test__str__(self, nucleotide_good):
        assert str(nucleotide_good) == 'A'


# class Test_letter_code_property  # tested in Test__Init__
# class Test_description_property
# class Test_degenerate_property
# class Test_supported_property
