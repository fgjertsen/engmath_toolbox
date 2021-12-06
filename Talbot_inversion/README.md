# Numerical inversion of Laplace transforms
Functions for inverse Laplace transforms, using  
1. Euler Inversion  
2. Talbot Inversion  

### List of example functions
0. Ramp function, type 1
1. Step function
2. Ramp function, type 2
3. Damped sine

## Usage
Build Euler Inversion function with `gcc -c euler_inversion.c`  
Build Talbot Inversion function with `gcc -c talbot_inversion.c`  

Build examples with `gcc -o Talbot ilt_examples.c talbot_inversion.o euler_inversion.o`  
Run with `Talbot i` where `i` is the index of the function to test.  
`Examples` will use `i = 3` as default.

Can be used directly with `Try_Talbot_Inversion.bat`.  

## Credits
**Main inspiration**: MATLAB code by Tucker McClure  
**Source**: Abate and Whitt (2006)