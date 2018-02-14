[![Build Status](https://travis-ci.org/komahanb/iterative-algebra.svg?branch=master)](https://travis-ci.org/komahanb/iterative-algebra) [![Coverage Status](https://coveralls.io/repos/github/komahanb/iterative-algebra/badge.svg?branch=master)](https://coveralls.io/github/komahanb/iterative-algebra?branch=master)

# Iterative Methods -- A Fortran-Python library for linear and nonlinear algebra

# How to install library?
sh ./build.sh

# How to run programs/examples?
- cd test
- cd specific_example
- make
- ./run

Directories:
------------
1. test/ contains the specific post-processing/type-setting done for coursework. Go to each sub-directory and 'make' and run the executable
2. src/ contains the Fortran source files. Run 'make; make install;'
3. python/ contains the python extension modules

Requirements:
------------
1. Modern Fortran compiler
2. Python 
3. f2py wrapper
