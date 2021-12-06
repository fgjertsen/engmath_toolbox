
/***************************************************************************************
* Header file for Talbot Inversion 													   *
***************************************************************************************/

#ifndef TALBOT_INVERSION_H
#define TALBOT_INVERSION_H

#include "ilt_common.h"

// Helper function: Cotangent
#define cot(x)	(cos(x)/sin(x))

// Main function prototype
double Talbot_inversion(					// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	);

#endif