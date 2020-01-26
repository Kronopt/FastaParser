@echo off

if "%1" == "" goto help
goto %~1

:help
echo.
echo install-dependencies   installs dependencies (includes development dependencies)
echo test                   runs tests
echo lint                   runs linter
echo coverage               runs test coverage
echo clean-tests            removes temp test files and folders
echo clean-coverage         removes coverage files
goto:eof

:install-dependencies
python -m pip install -r requirements-dev.txt
goto:eof

:test
python -m pytest tests/ -vv
goto:eof

:lint
python -m pylint fastaparser
goto:eof

:coverage
python -m coverage run --source fastaparser -m pytest tests/ -q
python -m coverage report -m
goto:eof

:clean-tests
rmdir /s /q .pytest_cache
goto:eof

:clean-coverage
python -m coverage erase
goto:eof
