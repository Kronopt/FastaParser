.PHONY: help install-dependencies test lint coverage build clean clean-pyc clean-tests clean-coverage clean-build

help:
	@echo ""
	@echo "install-dependencies     installs dependencies (includes development dependencies)"
	@echo "test                     runs tests"
	@echo "lint                     runs linter"
	@echo "coverage                 runs test coverage"
	@echo "build                    builds python package (sdist)"
	@echo ""
	@echo "clean                    runs all cleaning functions"
	@echo "clean-pyc                removes python file artifacts"
	@echo "clean-tests              removes temp test files and folders"
	@echo "clean-coverage           removes coverage files"
	@echo "clean-build              removes packaging artifacts"

install-dependencies:
	python -m pip install -r requirements-dev.txt

test:
	python -m pytest tests/ -vv

lint:
	python -m pylint fastaparser setup.py

coverage:
	python -m coverage run --source fastaparser -m pytest tests/ -q
	python -m coverage report -m

build: clean-pyc clean-build
	python setup.py sdist bdist_wheel

clean: clean-pyc clean-tests clean-coverage clean-build

clean-pyc:
	find . -name '__pycache__' -exec rm -rf {} +
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

clean-tests:
	rm -rf .pytest_cache/

clean-coverage:
	python -m coverage erase

clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf fastaparser.egg-info/
