
/***************************************************************************************
*	Inverse Laplace transforms														   *
*	Inspired by MATLAB code by Tucker McClure                                          *
*	Methods: 																		   *
*		1) Talbot's method, according to Abate and Whitt (2006)                        *
*		2) Euler's method															   *
*                                                                                      *
*	Fredrik Gjertsen 	<fredrik.gjertsen@gmail.com>  		                           *
***************************************************************************************/

#ifndef NILT_H
#define NILT_H

#include <stdlib.h>		/* Standard Library functions			*/
#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standard Library of Complex Numbers 	*/
#include <math.h>		/* Standard Library of Math				*/

#define M		64		// Algorithm parameter

// Helper function: Cotangent
#define cot(x)	(cos(x)/sin(x))

// Main function prototype: Euler inversion
double Euler_inversion(						// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	);

// Main function prototype: Talbot inversion
double Talbot_inversion(					// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	);

#endif