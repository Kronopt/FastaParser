.PHONY test

help:
    @echo "test               runs tests"
    @echo "clean-tests        removes temp test files and folders"

test:
	python -m pytest -v

clean-tests:
    rm -rf .pytest_cache/
