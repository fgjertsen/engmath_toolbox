# Numerical Inversion of Laplace transforms (__NILT__)
Functions for inverse Laplace transforms, using  
1. Euler Inversion  
2. Talbot Inversion  

### List of example functions
0. Ramp function, type 1
1. Step function
2. Ramp function, type 2
3. Damped sine

## Usage
Build Numerically Inversion functions with `gcc -c nilt.c`  

Build examples with `gcc -o NILT nilt_examples.c nilt.o`  
Run with `NILT i` where `i` is the index of the function to test.  
`NILT` will use `i = 3` as default.

Can be used directly with `Try_Talbot_Inversion.bat`.  

## Credits
**Main inspiration**: MATLAB code by Tucker McClure  
**Source**: Abate and Whitt (2006)

### Disclaimer
I did not invent these methods. I only built my own C-implementation of the methods.