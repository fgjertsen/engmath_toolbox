
/***************************************************************************************
* Header file for Euler Inversion													   *
***************************************************************************************/

#ifndef EULER_INVERSION_H
#define EULER_INVERSION_H

#include "ilt_common.h"

// Dimension for algorithm
static const int m = M/2;

// Main function prototype
double Euler_inversion(						// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	);

#endif