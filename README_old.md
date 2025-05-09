alphaMELTS build
================

Makefile_old was used for building the alphamelts executable and libalphamelts DLL before the C, Perl, and MATLAB/Python code were brought together into a single monolithic and open source repository. It is provided for reference so that required compiler flags etc. can be deduced until a poetry / cmake workflow is available.

Makefile_old is a recursive makefile. The .c source files in 'lib' were divided between the libphmelts and libalphamelts directories. The source code of the alphamelts executable (now in 'cli') were in a directory named 'alphamelts'.

Perl scripts in 'app' and MATLAB and Python source code were in separate repositories. The last version of these are at:
- https://magmasource.caltech.edu/gitlist/alphaMELTS2.git/
- https://magmasource.caltech.edu/gitlist/MELTS_Matlab.git/