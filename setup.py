#!python
# coding: utf-8

"""
FastaParser setup configuration.
"""

import fastaparser


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md', encoding='utf-8') as readme:
    DESCRIPTION = readme.read()
    readme.seek(0)
    next(readme)
    DESCRIPTION_SHORT = next(readme).rstrip()  # 2nd line of README.md is the short description

setup(
    name='fastaparser',
    version=fastaparser.__version__,
    description=DESCRIPTION_SHORT,
    long_description=DESCRIPTION,
    long_description_content_type='text/markdown',
    license=fastaparser.__license__,
    author=fastaparser.__author__,
    url='https://github.com/Kronopt/FastaParser',
    project_urls={'Documentation': 'https://fastaparser.readthedocs.io/en/latest/'},
    packages=['fastaparser'],
    package_dir={'fastaparser': 'fastaparser'},
    keywords='fasta parser',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],
    test_suite='tests',
    tests_require=['pytest'],
)
