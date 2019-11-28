#!python
# coding: utf-8

"""
Pytest configuration
"""


import pytest


@pytest.fixture(scope='session')
def unknown_characters():
    return 'O', '»', '%', 'º', '?', 'غ', '\n'
