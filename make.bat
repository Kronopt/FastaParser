@echo off

if "%1" == "" goto help
goto %~1

:help
echo.
echo test               runs tests
echo lint               run linter
echo clean-tests        removes temp test files and folders
goto:eof

:test
python -m pytest
goto:eof

:lint
python -m pylint FastaParser
goto:eof

:clean-tests
rmdir /s /q .pytest_cache
goto:eof
