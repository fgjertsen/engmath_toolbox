@echo off
REM Build Talbot Inversion functions
gcc -c talbot_inversion.c 
REM Build Euler Inversion functions
gcc -c euler_inversion.c

REM Build Comparison examples
gcc -o Talbot ilt_examples.c talbot_inversion.o euler_inversion.o
del *.o

REM Run examples
Talbot 3
pause