# Numerical inversion of Laplace transforms
Functions for inverse Laplace transforms, using
	1. Euler Inversion
	2. Talbot Inversion

### List of functions
0. Ramp function, type 1
1. Step function
2. Ramp function, type 2
3. Damped sine

## Usage
Build with `gcc talbot_inversion.c -o Talbot`.  
Run with `Talbot i` where `i` is the index of the function to test.  
`Talbot` will used `i=3` as default.  
Can be used directly with `Try_Talbot_Inversion.bat`.

## Credit
**Main inspiration**: MATLAB code by Tucker McClure
**Source**: Abate and Whitt (2006)