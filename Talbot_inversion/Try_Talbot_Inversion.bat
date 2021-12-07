@echo off
REM Build Numerical Inversion functions
gcc -c nilt.c 

REM Build Comparison examples
gcc -o NILT nilt_examples.c nilt.o
del *.o

REM Run examples
NILT 3
pause