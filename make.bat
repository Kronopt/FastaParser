@echo off

if "%1" == "" goto help
goto %~1

:help
echo.
echo install-dependencies   installs dependencies (includes development dependencies)
echo test                   runs tests
echo lint                   runs linter
echo clean-tests            removes temp test files and folders
goto:eof

:install-dependencies
python -m pip install -r requirements-dev.txt
goto:eof

:test
python -m pytest
goto:eof

:lint
python -m pylint fastaparser
goto:eof

:clean-tests
rmdir /s /q .pytest_cache
goto:eof
