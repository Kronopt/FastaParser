#!python
# coding: utf-8

"""
Tests for fastaparser.FastaSequence class.
"""


import pytest
from fastaparser import (
    FastaSequence,
    LetterCode,
    NUCLEOTIDE_LETTER_CODES_GOOD,
    AMINOACID_LETTER_CODES_GOOD,
    NUCLEOTIDE_LETTER_CODES_DEGENERATE,
    AMINOACID_LETTER_CODES_DEGENERATE,
    NUCLEOTIDE_LETTER_CODES_COMPLEMENT,
    AMINOACIDS_NOT_IN_NUCLEOTIDES,
)


##########
# Fixtures
##########


@pytest.fixture()
def nucleotide_good():
    return (
        FastaSequence(
            "".join(NUCLEOTIDE_LETTER_CODES_GOOD), sequence_type="nucleotide"
        ),
        [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD],
    )


@pytest.fixture()
def nucleotide_degenerate():
    return (
        FastaSequence(
            "".join(NUCLEOTIDE_LETTER_CODES_DEGENERATE), sequence_type="nucleotide"
        ),
        [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_DEGENERATE],
    )


@pytest.fixture()
def nucleotide_good_complement():
    return [
        LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT[letter_code])
        for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD
    ]


@pytest.fixture()
def aminoacid_good():
    return (
        FastaSequence("".join(AMINOACID_LETTER_CODES_GOOD), sequence_type="aminoacid"),
        [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_GOOD],
    )


@pytest.fixture()
def aminoacid_degenerate():
    return (
        FastaSequence(
            "".join(AMINOACID_LETTER_CODES_DEGENERATE), sequence_type="aminoacid"
        ),
        [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_DEGENERATE],
    )


@pytest.fixture()
def sequence_type_none():
    return (
        FastaSequence(
            "".join(NUCLEOTIDE_LETTER_CODES_GOOD) + "".join(AMINOACID_LETTER_CODES_GOOD)
        ),
        [LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD]
        + [LetterCode(letter_code) for letter_code in AMINOACID_LETTER_CODES_GOOD],
    )


@pytest.fixture()
def sequence_type_none_complement():
    return [
        LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT[letter_code])
        for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD
    ] + [
        LetterCode(NUCLEOTIDE_LETTER_CODES_COMPLEMENT.get(letter_code, letter_code))
        for letter_code in AMINOACID_LETTER_CODES_GOOD
    ]


@pytest.fixture()
def letter_codes_unknown(unknown_characters):
    return (
        FastaSequence("".join(unknown_characters)),
        [LetterCode(letter_code) for letter_code in unknown_characters],
    )


@pytest.fixture()
def actg_letter_code_list():
    return [LetterCode(letter_code) for letter_code in "ACTG"]


@pytest.fixture()
def actgnu_letter_code_counts():
    return {"A": 1, "C": 1, "G": 1, "T": 1, "N": 1, "U": 1}


#######
# Tests
#######


class Test__init__:
    def test_sequence_nucleotide_good(self, nucleotide_good):
        fasta_sequence, correct_sequence = nucleotide_good
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "nucleotide"
        assert fasta_sequence.inferred_type is False

    def test_sequence_nucleotide_degenerate(self, nucleotide_degenerate):
        fasta_sequence, correct_sequence = nucleotide_degenerate
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "nucleotide"
        assert fasta_sequence.inferred_type is False

    def test_sequence_aminoacid_good(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "aminoacid"
        assert fasta_sequence.inferred_type is False

    def test_sequence_aminoacid_degenerate(self, aminoacid_degenerate):
        fasta_sequence, correct_sequence = aminoacid_degenerate
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "aminoacid"
        assert fasta_sequence.inferred_type is False

    def test_sequence_unknown_characters(self, letter_codes_unknown):
        fasta_sequence, correct_sequence = letter_codes_unknown
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_sequence_not_str(self):
        with pytest.raises(TypeError):
            FastaSequence(1)
        with pytest.raises(TypeError):
            FastaSequence([])
        with pytest.raises(TypeError):
            FastaSequence(FastaSequence("actg"))

    def test_sequence_empty_str(self):
        with pytest.raises(TypeError):
            FastaSequence("")

    def test_sequence_lower_case(self):
        fasta_sequence = FastaSequence("".join(NUCLEOTIDE_LETTER_CODES_GOOD).lower())
        correct_sequence = [
            LetterCode(letter_code) for letter_code in NUCLEOTIDE_LETTER_CODES_GOOD
        ]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_id_good(self, actg_letter_code_list):
        id_ = "correct_id|some_other_id"
        # include '>'
        fasta_sequence = FastaSequence("ACTG", id_=">" + id_)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # exclude '>'
        fasta_sequence = FastaSequence("ACTG", id_=id_)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_id_with_spaces_and_newlines(self, actg_letter_code_list):
        id_ = "   correct id|some other id\n more id\r"
        fasta_sequence = FastaSequence("ACTG", id_=id_)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == "correct_id|some_other_id_more_id"
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_id_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG", id_=1)
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG", id_=[])

    # def test_id_empty_str (already tested)

    def test_description_good(self, actg_letter_code_list):
        description = "some correct description with numbers 1234567890#"
        fasta_sequence = FastaSequence("ACTG", description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_description_with_newlines(self, actg_letter_code_list):
        description = "   some correct description \nwith numbers \r1234567890#"
        fasta_sequence = FastaSequence("ACTG", description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert (
            fasta_sequence.description
            == "some correct description with numbers 1234567890#"
        )
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_description_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG", description=1)
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG", description=[])

    # def test_description_empty_str (already tested)

    def test_id_and_description_good(self, actg_letter_code_list):
        id_ = "correct_id|some_other_id"
        description = "some correct description with numbers 1234567890#"
        # include '>'
        fasta_sequence = FastaSequence("ACTG", id_=">" + id_, description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # exclude '>'
        fasta_sequence = FastaSequence("ACTG", id_=id_, description=description)
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    def test_sequence_type_nucleotide(self, unknown_characters, letter_codes_unknown):
        # nucleotide_good already tested
        fasta_sequence = FastaSequence(
            "".join(unknown_characters), sequence_type="nucleotide"
        )
        assert fasta_sequence.sequence == letter_codes_unknown[1]
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "nucleotide"
        assert fasta_sequence.inferred_type is False

    def test_sequence_type_aminoacid(self, unknown_characters, letter_codes_unknown):
        # aminoacid_good already tested
        fasta_sequence = FastaSequence(
            "".join(unknown_characters), sequence_type="aminoacid"
        )
        assert fasta_sequence.sequence == letter_codes_unknown[1]
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "aminoacid"
        assert fasta_sequence.inferred_type is False

    # def test_sequence_type_none (already tested)

    def test_sequence_type_wrong_type(self):
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type=1)
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type=[])

    def test_sequence_type_wrong_str(self):
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type="")
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type="wrong string")
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type="aminoacids")
        with pytest.raises(TypeError):
            FastaSequence("ACTG", sequence_type="nucleotides")

    def test_infer_type_true(self, nucleotide_good, aminoacid_good):
        # nucleotide sequence
        fasta_sequence = FastaSequence(
            "".join(NUCLEOTIDE_LETTER_CODES_GOOD), infer_type=True
        )
        correct_sequence = nucleotide_good[1]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # aminoacid sequence
        fasta_sequence = FastaSequence(
            "".join(AMINOACID_LETTER_CODES_GOOD), infer_type=True
        )
        correct_sequence = aminoacid_good[1]
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "aminoacid"
        assert fasta_sequence.inferred_type is True

    # def test_infer_type_false (already tested)

    def test_infer_type_not_bool(self):
        with pytest.raises(TypeError):
            FastaSequence("ACTG", infer_type=1)
        with pytest.raises(TypeError):
            FastaSequence("ACTG", infer_type=[])
        with pytest.raises(TypeError):
            FastaSequence("ACTG", infer_type="")

    def test_counts_good(self, nucleotide_good, aminoacid_good):
        assert nucleotide_good[0]._counts == dict(
            zip(NUCLEOTIDE_LETTER_CODES_GOOD, [1] * len(NUCLEOTIDE_LETTER_CODES_GOOD))
        )
        assert aminoacid_good[0]._counts == dict(
            zip(AMINOACID_LETTER_CODES_GOOD, [1] * len(AMINOACID_LETTER_CODES_GOOD))
        )

    def test_counts_same_character(self):
        sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        fasta_sequence = FastaSequence(sequence)
        assert fasta_sequence._counts == {"A": len(sequence)}

    def test_counts_letter_codes_unknown(
        self, letter_codes_unknown, unknown_characters
    ):
        assert letter_codes_unknown[0]._counts == dict(
            zip(unknown_characters, [1] * len(unknown_characters))
        )

    def test_gc(self, nucleotide_good):
        assert nucleotide_good[0]._gc is None

    def test_at(self, nucleotide_good):
        assert nucleotide_good[0]._at is None

    def test_current_iterator(self, nucleotide_good):
        assert nucleotide_good[0]._current_iterator is None


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
            FastaSequence.from_fastasequence("ACTG")
        with pytest.raises(TypeError):
            FastaSequence.from_fastasequence(1)


class Test_id_property:
    def test_set_good(self, actg_letter_code_list):
        id_ = "correct_id|some_other_id"
        # include '>'
        fasta_sequence = FastaSequence("ACTG")
        fasta_sequence.id = ">" + id_
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False
        # exclude '>'
        fasta_sequence = FastaSequence("ACTG")
        fasta_sequence.id = id_
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_good(already tested in Test__Init__)

    def test_set_with_spaces_and_newlines(self, actg_letter_code_list):
        id_ = "   correct id|some other id\n more id\r"
        fasta_sequence = FastaSequence("ACTG")
        fasta_sequence.id = id_
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == "correct_id|some_other_id_more_id"
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_with_spaces_and_newlines(already tested in Test__Init__)

    def test_set_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG")
            fasta_sequence.id = 1
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG")
            fasta_sequence.id = []

    # def test_get_not_str(already tested in Test__Init__)

    def test_delete(self, actg_letter_code_list):
        fasta_sequence = FastaSequence("ACTG", id_="test_id")
        del fasta_sequence.id
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False


class Test_description_property:
    def test_set_good(self, actg_letter_code_list):
        description = "some correct description with numbers 1234567890#"
        fasta_sequence = FastaSequence("ACTG")
        fasta_sequence.description = description
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == description
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_good(already tested in Test__Init__)

    def test_set_with_newlines(self, actg_letter_code_list):
        description = "   some correct description \nwith numbers \r1234567890#"
        fasta_sequence = FastaSequence("ACTG")
        fasta_sequence.description = description
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert (
            fasta_sequence.description
            == "some correct description with numbers 1234567890#"
        )
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_with_newlines(already tested in Test__Init__)

    def test_set_not_str(self):
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG")
            fasta_sequence.description = 1
        with pytest.raises(TypeError):
            fasta_sequence = FastaSequence("ACTG")
            fasta_sequence.description = []

    # def test_get_not_str(already tested in Test__Init__)

    def test_delete(self, actg_letter_code_list):
        fasta_sequence = FastaSequence("ACTG", description="test_description")
        del fasta_sequence.description
        assert fasta_sequence.sequence == actg_letter_code_list
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False


class Test_sequence_type_property:
    def test_set_nucleotide(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        fasta_sequence.sequence_type = "nucleotide"
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "nucleotide"
        assert fasta_sequence.inferred_type is False

    # def test_get_nucleotide (already tested in Test__Init__)

    def test_set_aminoacid(self, nucleotide_good):
        fasta_sequence, correct_sequence = nucleotide_good
        fasta_sequence.sequence_type = "aminoacid"
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type == "aminoacid"
        assert fasta_sequence.inferred_type is False

    # def test_get_aminoacid (already tested in Test__Init__)

    def test_set_none(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        fasta_sequence.sequence_type = None
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False

    # def test_get_none (already tested in Test__Init__)

    def test_set_wrong_str(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = "something"
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = ""

    # def test_get_wrong_str (already tested in Test__Init__)

    def test_set_wrong_type(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = 1
        with pytest.raises(TypeError):
            fasta_sequence.sequence_type = []

    # def test_get_wrong_type (already tested in Test__Init__)

    def test_delete(self, aminoacid_good):
        fasta_sequence, correct_sequence = aminoacid_good
        del fasta_sequence.sequence_type
        assert fasta_sequence.sequence == correct_sequence
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == ""
        assert fasta_sequence.sequence_type is None
        assert fasta_sequence.inferred_type is False


class Test_complement:
    def test_sequence_type_aminoacid(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.complement()

    def test_sequence_type_nucleotide_reverse_False(
        self, nucleotide_good, nucleotide_good_complement
    ):
        fasta_sequence = nucleotide_good[0]
        complement = fasta_sequence.complement()
        assert complement.sequence == nucleotide_good_complement
        assert complement.id == fasta_sequence.id
        assert complement.description == "[COMPLEMENT]"
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_sequence_type_nucleotide_reverse_true(
        self, nucleotide_good, nucleotide_good_complement
    ):
        fasta_sequence = nucleotide_good[0]
        complement = fasta_sequence.complement(True)
        assert complement.sequence == nucleotide_good_complement[::-1]
        assert complement.id == fasta_sequence.id
        assert complement.description == "[REVERSE COMPLEMENT]"
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_reverse_not_bool(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.complement(1)
        with pytest.raises(TypeError):
            fasta_sequence.complement([])
        with pytest.raises(TypeError):
            fasta_sequence.complement("")

    def test_sequence_type_none(
        self, sequence_type_none, sequence_type_none_complement
    ):
        fasta_sequence = sequence_type_none[0]
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == sequence_type_none_complement
        assert complement.id == fasta_sequence.id
        assert complement.description == "[COMPLEMENT]"
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type
        # aminoacid letter codes that are not also a nucleotide
        fasta_sequence = FastaSequence("".join(AMINOACIDS_NOT_IN_NUCLEOTIDES))
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == [
            LetterCode(letter_code) for letter_code in AMINOACIDS_NOT_IN_NUCLEOTIDES
        ]
        assert complement.id == fasta_sequence.id
        assert complement.description == "[COMPLEMENT]"
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type

    def test_letter_code_unknown(self, letter_codes_unknown):
        fasta_sequence, correct_sequence = letter_codes_unknown
        with pytest.warns(UserWarning):
            complement = fasta_sequence.complement()
        assert complement.sequence == correct_sequence
        assert complement.id == fasta_sequence.id
        assert complement.description == "[COMPLEMENT]"
        assert complement.sequence_type == fasta_sequence.sequence_type
        assert complement.inferred_type == fasta_sequence.inferred_type


class Test_gc_content:
    def test_sequence_letter_codes_gcs(self):
        for letter_code in "GCS":
            # single letter code sequence
            fasta_sequence = FastaSequence(letter_code, sequence_type="nucleotide")
            assert fasta_sequence._gc is None
            assert fasta_sequence.gc_content() == 1
            assert fasta_sequence._gc == 1
            assert fasta_sequence.gc_content(True) == 100
            assert fasta_sequence._gc == 1
            # multiple letter code sequence (only one G/C/S)
            fasta_sequence = FastaSequence(
                letter_code + "ATUATUATU", sequence_type="nucleotide"
            )
            assert fasta_sequence._gc is None
            assert fasta_sequence.gc_content() == 0.1
            assert fasta_sequence._gc == 1
            assert fasta_sequence.gc_content(True) == 10
            assert fasta_sequence._gc == 1

    def test_sequence_letter_codes_not_gcs(self):
        # single letter code sequence
        for letter_code in "ATU":
            fasta_sequence = FastaSequence(letter_code, sequence_type="nucleotide")
            assert fasta_sequence._gc is None
            assert fasta_sequence.gc_content() == 0
            assert fasta_sequence._gc == 0
            assert fasta_sequence.gc_content(True) == 0
            assert fasta_sequence._gc == 0
        # multiple letter code sequence
        fasta_sequence = FastaSequence("ATUATUATU", sequence_type="nucleotide")
        assert fasta_sequence._gc is None
        assert fasta_sequence.gc_content() == 0
        assert fasta_sequence._gc == 0
        assert fasta_sequence.gc_content(True) == 0
        assert fasta_sequence._gc == 0

    # def test_sequence_type_nucleotide (already tested)

    def test_sequence_type_aminoacid(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.gc_content()

    def test_sequence_type_none(self, sequence_type_none):
        fasta_sequence = sequence_type_none[0]
        with pytest.warns(UserWarning):
            fasta_sequence.gc_content()

    def test_gc_content_already_set(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        fasta_sequence._gc = 5  # would have been 2, by default
        # like this it shows the usage of the _gc variable, therefore no new gc_content is computed
        assert fasta_sequence.gc_content() == 5 / len(fasta_sequence.sequence)
        assert fasta_sequence._gc == 5
        # reset _gc
        fasta_sequence._gc = None
        assert fasta_sequence.gc_content() == 2 / len(fasta_sequence.sequence)
        assert fasta_sequence._gc == 2

    # def test_as_percentage_True (already tested)
    # def test_as_percentage_False (already tested)

    def test_as_percentage_not_bool(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.gc_content(1)
        with pytest.raises(TypeError):
            fasta_sequence.gc_content("")
        with pytest.raises(TypeError):
            fasta_sequence.gc_content([])


class Test_at_gc_ratio:
    def test_sequence_letter_codes_atw_gcs(self):
        for letter_code in "ATW":
            # single letter code sequence
            fasta_sequence = FastaSequence(letter_code, sequence_type="nucleotide")
            assert fasta_sequence._at is None
            assert fasta_sequence._gc is None
            assert fasta_sequence.at_gc_ratio() == 0
            assert fasta_sequence._at == 1
            assert fasta_sequence._gc == 0
            # multiple letter code sequence (only one A/T/W)
            fasta_sequence = FastaSequence(
                letter_code + "UNUNUNUNU", sequence_type="nucleotide"
            )
            assert fasta_sequence._at is None
            assert fasta_sequence._gc is None
            assert fasta_sequence.at_gc_ratio() == 0
            assert fasta_sequence._at == 1
            assert fasta_sequence._gc == 0
        for letter_code in "GCS":
            # single letter code sequence
            fasta_sequence = FastaSequence(letter_code, sequence_type="nucleotide")
            assert fasta_sequence._at is None
            assert fasta_sequence._gc is None
            assert fasta_sequence.at_gc_ratio() == 0
            assert fasta_sequence._at == 0
            assert fasta_sequence._gc == 1
            # multiple letter code sequence (only one G/C/S)
            fasta_sequence = FastaSequence(
                letter_code + "UNUNUNUNU", sequence_type="nucleotide"
            )
            assert fasta_sequence._at is None
            assert fasta_sequence._gc is None
            assert fasta_sequence.at_gc_ratio() == 0
            assert fasta_sequence._at == 0
            assert fasta_sequence._gc == 1
        # proportion 1:1
        fasta_sequence = FastaSequence("ATWGCS", sequence_type="nucleotide")
        assert fasta_sequence._at is None
        assert fasta_sequence._gc is None
        assert fasta_sequence.at_gc_ratio() == 1
        assert fasta_sequence._at == 3
        assert fasta_sequence._gc == 3
        # proportion 3:7
        fasta_sequence = FastaSequence("ATWGCSGCSG", sequence_type="nucleotide")
        assert fasta_sequence._at is None
        assert fasta_sequence._gc is None
        assert fasta_sequence.at_gc_ratio() == 3 / 7
        assert fasta_sequence._at == 3
        assert fasta_sequence._gc == 7
        # proportion 7:3
        fasta_sequence = FastaSequence("GCSATWATWA", sequence_type="nucleotide")
        assert fasta_sequence._at is None
        assert fasta_sequence._gc is None
        assert fasta_sequence.at_gc_ratio() == 7 / 3
        assert fasta_sequence._at == 7
        assert fasta_sequence._gc == 3

    def test_sequence_letter_codes_not_atw_gcs(self):
        # single letter code sequence
        for letter_code in "UN":
            fasta_sequence = FastaSequence(letter_code, sequence_type="nucleotide")
            assert fasta_sequence._at is None
            assert fasta_sequence._gc is None
            assert fasta_sequence.at_gc_ratio() == 0
            assert fasta_sequence._at == 0
            assert fasta_sequence._gc == 0
        # multiple letter code sequence
        fasta_sequence = FastaSequence("UNUNUNUNUN", sequence_type="nucleotide")
        assert fasta_sequence._at is None
        assert fasta_sequence._gc is None
        assert fasta_sequence.at_gc_ratio() == 0
        assert fasta_sequence._at == 0
        assert fasta_sequence._gc == 0

    # def test_sequence_type_nucleotide (already tested)

    def test_sequence_type_aminoacid(self, aminoacid_good):
        fasta_sequence = aminoacid_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.at_gc_ratio()

    def test_sequence_type_none(self, sequence_type_none):
        fasta_sequence = sequence_type_none[0]
        with pytest.warns(UserWarning):
            fasta_sequence.at_gc_ratio()

    def test_at_already_set(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        fasta_sequence._at = 9  # would have been 2, by default
        # like this it shows the usage of the _at variable, therefore no new at is computed
        assert fasta_sequence.at_gc_ratio() == 9 / 2
        assert fasta_sequence._at == 9
        assert fasta_sequence._gc == 2
        # reset _at and _gc
        fasta_sequence._at = None
        fasta_sequence._gc = None
        assert fasta_sequence.at_gc_ratio() == 1
        assert fasta_sequence._at == 2
        assert fasta_sequence._gc == 2

    def test_gc_already_set(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        fasta_sequence._gc = 3  # would have been 2, by default
        # like this it shows the usage of the _gc variable, therefore no new gc is computed
        assert fasta_sequence.at_gc_ratio() == 2 / 3
        assert fasta_sequence._at == 2
        assert fasta_sequence._gc == 3
        # reset _at and _gc
        fasta_sequence._at = None
        fasta_sequence._gc = None
        assert fasta_sequence.at_gc_ratio() == 1
        assert fasta_sequence._at == 2
        assert fasta_sequence._gc == 2

    def test_at_and_gc_already_set(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        fasta_sequence._at = 9  # would have been 2, by default
        fasta_sequence._gc = 3  # would have been 2, by default
        # like this it shows the usage of the _at and _gc variables, therefore no new at or gc is computed
        assert fasta_sequence.at_gc_ratio() == 9 / 3
        assert fasta_sequence._at == 9
        assert fasta_sequence._gc == 3
        # reset _at and _gc
        fasta_sequence._at = None
        fasta_sequence._gc = None
        assert fasta_sequence.at_gc_ratio() == 1
        assert fasta_sequence._at == 2
        assert fasta_sequence._gc == 2


class Test_count_letter_codes:
    def test_letter_codes_None(self, nucleotide_good, actgnu_letter_code_counts):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        assert fasta_sequence.count_letter_codes() == fasta_sequence._counts
        assert fasta_sequence.count_letter_codes() == actgnu_letter_code_counts

    def test_letter_codes_iterable(self, nucleotide_good, actgnu_letter_code_counts):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        # empty iterables
        assert fasta_sequence.count_letter_codes("") == actgnu_letter_code_counts
        assert fasta_sequence.count_letter_codes("") == fasta_sequence._counts
        assert fasta_sequence.count_letter_codes([]) == actgnu_letter_code_counts
        assert fasta_sequence.count_letter_codes([]) == fasta_sequence._counts
        assert fasta_sequence.count_letter_codes({}) == actgnu_letter_code_counts
        assert fasta_sequence.count_letter_codes({}) == fasta_sequence._counts
        assert fasta_sequence.count_letter_codes(()) == actgnu_letter_code_counts
        assert fasta_sequence.count_letter_codes(()) == fasta_sequence._counts
        # iterables
        del actgnu_letter_code_counts["U"]  # just to differ from the default
        assert fasta_sequence.count_letter_codes("ACGTN") == actgnu_letter_code_counts
        assert (
            fasta_sequence.count_letter_codes(["A", "C", "G", "T", "N"])
            == actgnu_letter_code_counts
        )
        assert (
            fasta_sequence.count_letter_codes(
                {"A": "blaablaaa", "C": 2112, "G": [], "T": "T", "N": ()}
            )
            == actgnu_letter_code_counts
        )
        assert (
            fasta_sequence.count_letter_codes(("A", "C", "G", "T", "N"))
            == actgnu_letter_code_counts
        )

    # def test_letter_codes_none (already tested)

    def test_letter_code_not_iterable(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.count_letter_codes(1111)
        with pytest.raises(TypeError):
            fasta_sequence.count_letter_codes(2.333)
        with pytest.raises(TypeError):
            fasta_sequence.count_letter_codes(True)

    # def test_letter_codes_in_sequence (already tested)

    def test_letter_codes_not_in_sequence(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]  # ACGTNU
        assert fasta_sequence.count_letter_codes("X") == {"X": 0}
        assert fasta_sequence.count_letter_codes("XYZ") == {"X": 0, "Y": 0, "Z": 0}
        assert fasta_sequence.count_letter_codes("ACXYZ") == {
            "A": 1,
            "C": 1,
            "X": 0,
            "Y": 0,
            "Z": 0,
        }


class Test_count_letter_codes_degenerate:
    def test_sequence_type_nucleotide(self, nucleotide_good, nucleotide_degenerate):
        fasta_sequence_good = nucleotide_good[0]
        fasta_sequence_degenerate = nucleotide_degenerate[0]
        assert fasta_sequence_good.count_letter_codes_degenerate() == {}
        assert fasta_sequence_degenerate.count_letter_codes_degenerate() == {
            "K": 1,
            "S": 1,
            "Y": 1,
            "M": 1,
            "W": 1,
            "R": 1,
            "B": 1,
            "D": 1,
            "H": 1,
            "V": 1,
            "-": 1,
        }

    def test_sequence_type_aminoacid(self, aminoacid_good, aminoacid_degenerate):
        fasta_sequence_good = aminoacid_good[0]
        fasta_sequence_degenerate = aminoacid_degenerate[0]
        assert fasta_sequence_good.count_letter_codes_degenerate() == {}
        assert fasta_sequence_degenerate.count_letter_codes_degenerate() == {"-": 1}

    def test_sequence_type_none(self, sequence_type_none):
        fasta_sequence = sequence_type_none[0]
        with pytest.raises(TypeError):
            fasta_sequence.count_letter_codes_degenerate()


class Test_formatted_definition_line:
    def test_id_and_description(self):
        id_ = "id123|id456"
        description = "this is a description"
        # with '>'
        fasta_sequence = FastaSequence("A", id_=">" + id_, description=description)
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.formatted_definition_line() == ">%s %s" % (
            id_,
            description,
        )
        # without '>'
        fasta_sequence = FastaSequence("A", id_=id_, description=description)
        assert fasta_sequence.id == id_
        assert fasta_sequence.description == description
        assert fasta_sequence.formatted_definition_line() == ">%s %s" % (
            id_,
            description,
        )

    def test_empty_id(self):
        description = "this is a description"
        fasta_sequence = FastaSequence("A", description=description)
        assert fasta_sequence.id == ""
        assert fasta_sequence.description == description
        assert fasta_sequence.formatted_definition_line() == "> " + description

    def test_empty_description(self):
        id_ = "id123|id456"
        fasta_sequence = FastaSequence("A", id_=id_)
        assert fasta_sequence.id == id_
        assert fasta_sequence.formatted_definition_line() == ">" + id_

    def test_empty_id_and_description(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence.formatted_definition_line() == ">"

    def test_only_one_line(self):
        id_ = "id123 \nab \rcd"
        description = "description with \n new lines \r123"
        fasta_sequence = FastaSequence("A", id_=id_, description=description)
        assert fasta_sequence.id == "id123_ab_cd"
        assert fasta_sequence.description == "description with new lines 123"
        assert (
            fasta_sequence.formatted_definition_line()
            == ">id123_ab_cd description with new lines 123"
        )


class Test_formatted_sequence:
    def test_good(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence.formatted_sequence() == "".join(
            NUCLEOTIDE_LETTER_CODES_GOOD
        )
        fasta_sequence_71A = FastaSequence("A" * 71)
        assert fasta_sequence_71A.formatted_sequence() == "A" * 70 + "\nA"

    def test_max_characters_per_line_custom_ranges(self):
        fasta_sequence = FastaSequence("A" * 10)
        assert (
            fasta_sequence.formatted_sequence(10)
            == "A" * 10
            == fasta_sequence.sequence_as_string()
        )
        assert (
            fasta_sequence.formatted_sequence(11)
            == "A" * 10
            == fasta_sequence.sequence_as_string()
        )
        assert fasta_sequence.formatted_sequence(4) == "%s\n%s\nAA" % ("A" * 4, "A" * 4)
        assert fasta_sequence.formatted_sequence(5) == "%s\n%s" % ("A" * 5, "A" * 5)

    def test_max_characters_per_line_zero_or_less(self, nucleotide_good):
        correct_format = "\n".join(NUCLEOTIDE_LETTER_CODES_GOOD)
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence.formatted_sequence(0) == correct_format
        assert fasta_sequence.formatted_sequence(-1) == correct_format

    def test_max_characters_per_line_9999999(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence.formatted_sequence(9999999) == "".join(
            NUCLEOTIDE_LETTER_CODES_GOOD
        )

    def test_max_characters_per_line_not_int(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence.formatted_sequence("")
        with pytest.raises(TypeError):
            fasta_sequence.formatted_sequence([])
        with pytest.raises(TypeError):
            fasta_sequence.formatted_sequence(LetterCode)


class Test_formatted_fasta:
    def test_good(self):
        fasta_sequence = FastaSequence(
            "A" * 150,
            id_="  id | 123 \n456",
            description="  this is a description\n 123 ",
        )
        assert fasta_sequence.formatted_fasta() == ">%s %s\n%s" % (
            "id_|_123_456",
            "this is a description 123",
            "%s\n%s\n%s" % ("A" * 70, "A" * 70, "A" * 10),
        )


class Test_sequence_as_string:
    def test_good(
        self, nucleotide_good, aminoacid_good, letter_codes_unknown, unknown_characters
    ):
        nucleotide_sequence = nucleotide_good[0]
        aminoacid_sequence = aminoacid_good[0]
        unknown_sequence = letter_codes_unknown[0]
        assert nucleotide_sequence.sequence_as_string() == "".join(
            NUCLEOTIDE_LETTER_CODES_GOOD
        )
        assert aminoacid_sequence.sequence_as_string() == "".join(
            AMINOACID_LETTER_CODES_GOOD
        )
        assert unknown_sequence.sequence_as_string() == "".join(unknown_characters)


class Test_reverse:
    def test_good(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        sequence = "".join(map(str, fasta_sequence.reverse()))
        assert sequence == fasta_sequence.sequence_as_string()[::-1]

    def test_is_iterable(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        try:
            iter(fasta_sequence.reverse())
        except TypeError:
            pytest.fail("FastaSequence.reverse() is not iterable.")


class Test_update_sequence_type:
    def test_update_sequence_type_update_letter_code_objects_not_bool(
        self, nucleotide_good
    ):
        fasta_sequence = nucleotide_good[0]
        with pytest.raises(TypeError):
            fasta_sequence._update_sequence_type(
                "aminoacid", update_letter_code_objects=" "
            )
        with pytest.raises(TypeError):
            fasta_sequence._update_sequence_type(
                "aminoacid", update_letter_code_objects=[]
            )

    # def test_update_sequence_type_good (already tested)
    # def test_update_sequence_type_wrong_type (already tested)


class Test__iter__:
    def test__iter__(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        iterated_sequence = [letter_code for letter_code in fasta_sequence.__iter__()]
        assert iterated_sequence == fasta_sequence._sequence
        assert iterated_sequence == fasta_sequence.sequence
        try:
            iter(fasta_sequence._current_iterator)
        except TypeError:
            pytest.fail(
                "When calling __iter__, FastaSequence._current_iterator is not iterable."
            )

    def test_iterate(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        iterated_sequence = [letter_code for letter_code in fasta_sequence]
        assert iterated_sequence == fasta_sequence._sequence
        assert iterated_sequence == fasta_sequence.sequence
        try:
            iter(fasta_sequence._current_iterator)
        except TypeError:
            pytest.fail(
                "When iterating, FastaSequence._current_iterator is not iterable."
            )

    def test_is_iterable(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        try:
            iter(fasta_sequence.__iter__())
        except TypeError:
            pytest.fail("FastaSequence.__iter__() is not iterable.")


class Test__reversed__:
    def test__reversed__(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        iterated_sequence = [
            letter_code for letter_code in fasta_sequence.__reversed__()
        ]
        assert iterated_sequence == fasta_sequence._sequence[::-1]
        assert iterated_sequence == fasta_sequence.sequence[::-1]
        try:
            iter(fasta_sequence._current_iterator)
        except TypeError:
            pytest.fail(
                "When calling __reversed__, FastaSequence._current_iterator is not iterable."
            )

    def test_iterate(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        iterated_sequence = [letter_code for letter_code in reversed(fasta_sequence)]
        assert iterated_sequence == fasta_sequence._sequence[::-1]
        assert iterated_sequence == fasta_sequence.sequence[::-1]
        try:
            iter(fasta_sequence._current_iterator)
        except TypeError:
            pytest.fail(
                "When iterating in reverse, FastaSequence._current_iterator is not iterable."
            )

    def test_is_iterable(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        try:
            iter(fasta_sequence.__reversed__())
        except TypeError:
            pytest.fail("FastaSequence.__reversed__() is not iterable.")


class Test__next__:
    def test_next_existing_iterator(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        iter(fasta_sequence)
        assert (
            fasta_sequence.__next__()
            == fasta_sequence.sequence[0]
            == nucleotide_good[1][0]
        )
        assert (
            next(fasta_sequence) == fasta_sequence.sequence[1] == nucleotide_good[1][1]
        )
        assert (
            next(fasta_sequence._current_iterator)
            == fasta_sequence.sequence[2]
            == nucleotide_good[1][2]
        )

    def test_next_non_existing_iterator(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert (
            fasta_sequence.__next__()
            == fasta_sequence.sequence[0]
            == nucleotide_good[1][0]
        )
        assert (
            next(fasta_sequence) == fasta_sequence.sequence[1] == nucleotide_good[1][1]
        )
        assert (
            next(fasta_sequence._current_iterator)
            == fasta_sequence.sequence[2]
            == nucleotide_good[1][2]
        )
        try:
            iter(fasta_sequence._current_iterator)
        except TypeError:
            pytest.fail(
                "When calling __next__, FastaSequence._current_iterator is not iterable."
            )


class Test__getitem__:
    def test_get_element_at_first_position(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        fasta_sequence_0 = fasta_sequence[0]
        assert fasta_sequence_0 == nucleotide_good[1][0] == fasta_sequence._sequence[0]

    def test_get_element_at_last_position(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        fasta_sequence_last = fasta_sequence[-1]
        assert (
            fasta_sequence_last
            == nucleotide_good[1][-1]
            == fasta_sequence._sequence[-1]
        )

    def test_get_slice_single(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]

        fasta_sequence_0_sliced = fasta_sequence[:1]
        assert (
            fasta_sequence_0_sliced._sequence
            == [nucleotide_good[1][0]]
            == fasta_sequence._sequence[:1]
        )
        assert fasta_sequence_0_sliced.description == "[SLICE OF ORIGINAL: None:1:None]"
        fasta_sequence_last_sliced = fasta_sequence[-1:]
        assert (
            fasta_sequence_last_sliced._sequence
            == [nucleotide_good[1][-1]]
            == fasta_sequence._sequence[-1:]
        )
        assert (
            fasta_sequence_last_sliced.description
            == "[SLICE OF ORIGINAL: -1:None:None]"
        )

        fasta_sequence_with_description_0_sliced = FastaSequence(
            "ACTG", description="existing description"
        )[:1]
        assert (
            fasta_sequence_with_description_0_sliced._sequence
            == ["A"]
            == fasta_sequence_with_description_0_sliced._sequence[:1]
        )
        assert fasta_sequence_with_description_0_sliced.description == (
            "existing description [SLICE OF ORIGINAL: " "None:1:None]"
        )
        fasta_sequence_description_last_sliced = FastaSequence(
            "ACTG", description="existing description"
        )[-1:]
        assert (
            fasta_sequence_description_last_sliced._sequence
            == ["G"]
            == fasta_sequence_description_last_sliced._sequence[-1:]
        )
        assert fasta_sequence_description_last_sliced.description == (
            "existing description [SLICE OF ORIGINAL: " "-1:None:None]"
        )

    def test_get_slice_multiple(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]

        fasta_sequence_sliced = fasta_sequence[2:4]
        assert (
            fasta_sequence_sliced._sequence
            == nucleotide_good[1][2:4]
            == fasta_sequence._sequence[2:4]
        )
        assert fasta_sequence_sliced.description == "[SLICE OF ORIGINAL: 2:4:None]"

    def test_get_slice_empty(self):
        with pytest.raises(TypeError):
            FastaSequence("ACTG")[:0]

    def test_get_item_wrong_type(self):
        with pytest.raises(TypeError):
            FastaSequence("ACTG").__getitem__("A")
        with pytest.raises(TypeError):
            FastaSequence("ACTG").__getitem__([])
        with pytest.raises(TypeError):
            FastaSequence("ACTG").__getitem__(LetterCode)


class Test__eq__:
    def test_other_fastasequence(self, nucleotide_good, aminoacid_good):
        fasta_sequence_nucleotide = nucleotide_good[0]
        fasta_sequence_aminoacid = aminoacid_good[0]
        fasta_sequence_nucleotide_equal = FastaSequence(
            "".join(NUCLEOTIDE_LETTER_CODES_GOOD)
        )
        assert fasta_sequence_nucleotide == fasta_sequence_nucleotide_equal
        assert fasta_sequence_nucleotide != fasta_sequence_aminoacid

    def test_other_str(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence == "".join(NUCLEOTIDE_LETTER_CODES_GOOD)
        assert fasta_sequence != "".join(NUCLEOTIDE_LETTER_CODES_GOOD)[::-1]
        assert fasta_sequence != "".join(NUCLEOTIDE_LETTER_CODES_GOOD)[:-1]
        assert fasta_sequence != ""

    def test_other_list(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence == [
            LetterCode(letter_code)
            for letter_code in "".join(NUCLEOTIDE_LETTER_CODES_GOOD)
        ]
        assert fasta_sequence != [
            LetterCode(letter_code)
            for letter_code in "".join(NUCLEOTIDE_LETTER_CODES_GOOD)[::-1]
        ]
        assert fasta_sequence != [
            LetterCode(letter_code)
            for letter_code in "".join(NUCLEOTIDE_LETTER_CODES_GOOD)[:-1]
        ]
        assert fasta_sequence != []

    def test_other_wrong_type(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert fasta_sequence != 1
        assert fasta_sequence != 0
        assert fasta_sequence != {
            letter_code: LetterCode(letter_code)
            for letter_code in "".join(NUCLEOTIDE_LETTER_CODES_GOOD)
        }
        assert fasta_sequence != {}
        assert fasta_sequence != LetterCode("A")


class Test__len__:
    def test__len__(self, nucleotide_good, aminoacid_good):
        fasta_sequence_nucleotide = nucleotide_good[0]
        fasta_sequence_aminoacid = aminoacid_good[0]
        assert len(fasta_sequence_nucleotide) == 6
        assert len(fasta_sequence_aminoacid) == 25


class Test__repr__:
    def test__repr__(self, nucleotide_good):
        fasta_sequence = nucleotide_good[0]
        assert repr(fasta_sequence) == "FastaSequence('%s')" % "".join(
            NUCLEOTIDE_LETTER_CODES_GOOD
        )


class Test__str__:
    def test__str__(self):
        fasta_sequence = FastaSequence("ACGTNU", id_="id", description="description")
        assert str(fasta_sequence) == ">id description\nACGTNU"


# tested in Test__Init__:
#   class Test_id_property
#   class Test_description_property
#   class Test_sequence_property
#   class Test_inferred_type
