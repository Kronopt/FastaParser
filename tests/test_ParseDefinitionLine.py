#!python
# coding: utf-8

"""
Tests for fastaparser.ParseDefinitionLine class.
"""


import pytest
from fastaparser import ParseDefinitionLine


##########
# Fixtures
##########


@pytest.fixture()
def definition_line_test_function():
    def definition_line_test(definition_line, final_id, final_description):
        id_, description = ParseDefinitionLine._parse_definition_line(definition_line)
        assert id_ == final_id
        assert description == final_description
    return definition_line_test


#######
# Tests
#######


class Test_parse_definition_line:
    def test_definition_line_id_and_description_good(self, definition_line_test_function):
        definition_line = '>ID123|moreID description and more text'
        definition_line_test_function(definition_line, 'ID123|moreID', 'description and more text')
        # without '>'
        definition_line_test_function(definition_line[1:], 'ID123|moreID', 'description and more text')
        # long space before description
        definition_line_test_function('id              description abc', 'id', 'description abc')

    def test_definition_line_id(self, definition_line_test_function):
        definition_line = '>ID123|secondID|otherID'
        definition_line_test_function(definition_line, 'ID123|secondID|otherID', '')
        # without '>'
        definition_line_test_function(definition_line[1:], 'ID123|secondID|otherID', '')
        # space after id
        definition_line_test_function(definition_line[1:] + '   ', 'ID123|secondID|otherID', '')

    def test_definition_line_empty(self, definition_line_test_function):
        definition_line_test_function('', '', '')
        # with '>'
        definition_line_test_function('>', '', '')

    def test_definition_line_wrong_type(self):
        with pytest.raises(TypeError):
            ParseDefinitionLine._parse_definition_line(1)
        with pytest.raises(TypeError):
            ParseDefinitionLine._parse_definition_line([])
