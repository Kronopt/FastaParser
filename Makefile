.PHONY: help install-dependencies test clean-tests

help:
	@echo "install-dependencies     installs dependencies (includes development dependencies)"
	@echo "test                     runs tests"
	@echo "lint                     runs linter"
	@echo "coverage                 runs test coverage"
	@echo "clean-tests              removes temp test files and folders"
	@echo "clean-coverage           removes coverage files"

install-dependencies:
	python -m pip install -r requirements-dev.txt

test:
	python -m pytest tests/ -vv

lint:
	python -m pylint fastaparser

coverage:
    python -m coverage run --source fastaparser -m pytest tests/ -q
    python -m coverage report -m

clean-tests:
	rm -rf .pytest_cache/

clean-coverage:
    rm -f .coverage
