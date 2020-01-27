@echo off

if "%1" == "" goto help
goto %~1

:help
echo.
echo install-dependencies   installs dependencies (includes development dependencies)
echo test                   runs tests
echo lint                   runs linter
echo coverage               runs test coverage
echo build                  builds python package (sdist)
echo build-test             tests build for errors and uploads to test.pypi.org
echo release                builds and uploads python package to pypi.org
echo.
echo clean                  runs all cleaning functions
echo clean-pyc              removes python file artifacts
echo clean-tests            removes temp test files and folders
echo clean-coverage         removes coverage files
echo clean-build            removes packaging artifacts
goto:eof

:install-dependencies
python -m pip install -r requirements-dev.txt
goto:eof

:test
python -m pytest tests/ -vv
goto:eof

:lint
python -m pylint fastaparser setup.py
goto:eof

:coverage
python -m coverage run --source fastaparser -m pytest tests/ -q
python -m coverage report -m
goto:eof

:build
call:clean-pyc
call:clean-build
python setup.py sdist bdist_wheel
goto:eof

:build-test
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
goto:eof

:release
call:build
twine upload dist/*
goto:eof

:clean
call:clean-pyc
call:clean-tests
call:clean-coverage
call:clean-build
goto:eof

:clean-pyc
rmdir /s /q fastaparser\__pycache__
rmdir /s /q tests\__pycache__
del /s *.pyc *.pyo *~
goto:eof

:clean-tests
rmdir /s /q .pytest_cache
goto:eof

:clean-coverage
python -m coverage erase
goto:eof

:clean-build
rmdir /s /q build
rmdir /s /q dist
rmdir /s /q fastaparser.egg-info
goto:eof
