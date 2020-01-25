.PHONY: help install-dependencies test clean-tests

help:
    @echo "install-dependencies     installs dependencies (includes development dependencies)"
    @echo "test                     runs tests"
    @echo "lint                     runs linter"
    @echo "clean-tests        removes temp test files and folders"

install-dependencies:
    python -m pip install -r requirements-dev.txt

test:
	python -m pytest

lint:
    python -m pylint FastaParser

clean-tests:
	rm -rf .pytest_cache/
