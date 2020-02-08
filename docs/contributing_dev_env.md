# Setting up Development Environment
To setup the development environment for FastaParser:

**1.** First `fork` FastaParser Github repository.

**2.** Then `clone` your fork into your machine. You can use [Github Desktop](https://desktop.github.com/) or git like so:
```sh
$ git clone git@github.com:your_name_here/FastaParser.git
```
Alternatively, you can simply download FastaParser from Github as a .zip file.

**3.** Although FastaParser targets Python versions 3.5-3.8, I recommend `Python 3.7` for development. As of the writing
of this documentation, mkdocs serve command won't work on python 3.8.

**4.** Create a virtual environment. Most IDE's have an option to create virtual environments.
Otherwise you can use virtualenv to create one.

Assuming you have virtualenv installed:
```sh
$ mkdir venv
$ virtualenv venv/FastaParser
```

to activate the enviroment:
```sh
$ cd venv/FastaParser/bin
$ source activate
```

**5.** Install the required packages for development. Run the following on the root of the FastaParser project:
```sh
$ make install-dependencies
```

**6.** You can now make all the changes you want.

**7.** When you're satisfied with your changes, push them to your fork of FastaParser. Again, you can use
[Github Desktop](https://desktop.github.com/) or git:
```sh
$ git add .
$ git commit -m "Your detailed description of your changes"
$ git push origin name-of-your-bugfix-or-feature
```

## Pull Requests
Before submitting your changes as a pull request:

**1.** Test your code by running the following on the root of the FastaParser project:
```sh
$ make test
$ make lint
$ make coverage
```
All tests should pass, linting should be 10/10 and coverage should be 100% or close.

If you added any documentation, check if it looks ok by running:
```sh
$ make docs-test
```
and opening the locally served documentation (you can close it by hitting Ctrl-C on the command line).

**2.** Don't forget to write docstrings, tests and documentation, if appropriate.
Try to follow the overall coding style of FastaParser.

**3.** Finally submit a pull request on the [FastaParser Github repository](https://github.com/Kronopt/FastaParser/pulls).
