.PHONY: help test clean-tests

help:
	@echo "test               runs tests"
	@echo "lint               run linter"
	@echo "clean-tests        removes temp test files and folders"

test:
	python -m pytest

lint:
    python -m pylint FastaParser

clean-tests:
	rm -rf .pytest_cache/
