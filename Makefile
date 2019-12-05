.PHONY: help test clean-tests

help:
    @echo "test               runs tests"
    @echo "clean-tests        removes temp test files and folders"

test:
	python -m pytest

clean-tests:
    rm -rf .pytest_cache/
