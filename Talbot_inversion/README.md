# Numerical Inversion of Laplace transforms (*NILT*)
Functions for inverse Laplace transforms, using  
1. Euler Inversion  
2. Talbot Inversion  

### List of example functions
0. Ramp function, type 1
1. Step function
2. Ramp function, type 2
3. Damped sine

## Usage
Build Numerical Inversion (of Laplace Transforms) functions with `gcc -c nilt.c`  

Build examples with `gcc -o NILT nilt_examples.c nilt.o`  
Run with `NILT i` where `i` is the index of the function to test.  
`NILT` will use `i = 3` as default.

Can be used directly with `Try_Talbot_Inversion.bat`.  

![I do not own rights to this picture](assets/P_S_Laplace.jpg?raw=true "Laplace")  
*(drawing acquired from Wikimedia Commons)*

## Credits
**Inspiration, implementation**: [MATLAB code by Tucker McClure](https://se.mathworks.com/matlabcentral/fileexchange/39035-numerical-inverse-laplace-transform)  
**Theory, methods**: [Abate and Whitt (2006)](http://www.columbia.edu/~ww2040/UnifiedDraft.pdf)  

### Disclaimer
Software "as is". No warranty.  
I did not invent these methods. I only made my own C-implementation of the methods.