.PHONY: help install-dependencies install-testingdependencies test lint coverage docs-test build build-test release clean clean-pyc clean-tests clean-coverage clean-build conda-install-dependencies conda-skeleton cconda-build-and-upload conda-clean-build

help:
	@echo ""
	@echo "install-dependencies         installs dependencies (includes development dependencies)"
	@echo "install-testing-dependencies installs dependencies for testing"
	@echo "test                         runs tests"
	@echo "lint                         runs linter"
	@echo "coverage                     runs test coverage"
	@echo "docs-test                    tests docs for build errors and serves them locally"
	@echo "build                        builds python package (sdist)"
	@echo "build-test                   tests build for errors and uploads to test.pypi.org"
	@echo "release                      builds and uploads python package to pypi.org"
	@echo ""
	@echo "clean                        runs all cleaning functions"
	@echo "clean-pyc                    removes python file artifacts"
	@echo "clean-tests                  removes temp test files and folders"
	@echo "clean-coverage               removes coverage files"
	@echo "clean-build                  removes packaging artifacts"
	@echo ""
	@echo "conda-install-dependencies   installs conda build dependencies"
	@echo "conda-skeleton               creates skeleton conda package recipe"
	@echo "conda-build-and-upload       builds conda package and uploads to anaconda cloud"
	@echo ""
	@echo "conda-clean-build            removes conda build artifacts"

install-dependencies:
	python -m pip install -r requirements-dev.txt

install-testing-dependencies:
	python -m pip install -r requirements-testing.txt

test:
	python -m pytest tests/ -vv

lint:
	python -m pylint fastaparser setup.py
	python -m black --check .

coverage:
	python -m coverage run --source fastaparser -m pytest tests/ -q
	python -m coverage report -m

docs-test:
	mkdocs serve -s -f .mkdocs.yml

build: clean-pyc clean-build
	python setup.py sdist bdist_wheel

build-test:
	twine check dist/*
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

release: build
	twine upload dist/*

clean: clean-pyc clean-tests clean-coverage clean-build

clean-pyc:
	rm -rf fastaparser/__pycache__ tests/__pycache__
	rm -f fastaparser/*.pyc tests/*.pyc
	rm -f fastaparser/*.pyo tests/*.pyo
	rm -f fastaparser/*~ tests/*~

clean-tests:
	rm -rf .pytest_cache/

clean-coverage:
	python -m coverage erase

clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf fastaparser.egg-info/

conda-install-dependencies:
	conda install conda-build anaconda-client

conda-skeleton:
	conda skeleton pypi fastaparser --output-dir .conda

conda-build-and-upload:
	conda config --set anaconda_upload yes
	conda build .conda --strict-verify --output-folder .build_conda

conda-clean-build:
	conda build purge
	rm -rf .build_conda/
