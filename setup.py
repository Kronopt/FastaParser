#!python
# coding: utf-8

"""
FastaParser setup configuration.
"""

from setuptools import setup
import fastaparser


with open("README.md", encoding="utf-8") as readme, open(
    "docs/history.md", encoding="utf-8"
) as history:
    DESCRIPTION = readme.read() + "\n#" + history.read()
DESCRIPTION_SHORT = "A Python FASTA file Parser and Writer."

setup(
    name="fastaparser",
    version=fastaparser.__version__,
    description=DESCRIPTION_SHORT,
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    license=fastaparser.__license__,
    author=fastaparser.__author__,
    url="https://github.com/Kronopt/FastaParser",
    project_urls={"Documentation": "https://fastaparser.readthedocs.io/en/latest/"},
    packages=["fastaparser"],
    package_dir={"fastaparser": "fastaparser"},
    keywords="fasta parser fasta-parser fasta-writer",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    test_suite="tests",
    tests_require=["pytest"],
)
