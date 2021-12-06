
/***************************************************************************************
*	Inverse Laplace transform, Talbot Inversion Method								   *
*	Inspired by MATLAB code by Tucker McClure                                          *
*   Method according to Abate and Whitt (2006)                                         *
*	Fredrik Gjertsen 	<fredrik.gjertsen@gmail.com>  		                           *
***************************************************************************************/

#include "talbot_inversion.h"

double Talbot_inversion(					// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	)
{
	double complex
		delta[M],
		gamma[M],
		gamma_fs[M];

	// Calculate delta
	delta[0] = 2.0*M/5.0;
	for (int k = 1; k < M; k++)
		delta[k]	= 2.0*M_PI/5.0 * k * ( cot(M_PI/M*k) + I );

	// Calculate gamma
	gamma[0] = 0.5;
	for (int k = 1; k < M; k++){
		gamma[k]	= M_PI/M*k*(1.0 + cot(M_PI/M*k)*cot(M_PI/M*k)) - cot(M_PI/M*k);
		gamma[k]   *= I;
		gamma[k]   += 1.0;
	}
	for (int k = 0; k < M; k++)
		gamma[k]   *= cexp(delta[k]);
	
	// Calculate inverse Laplace Transform
	double ilt = 0.0;
	for (int k = 0; k < M; k++){
		gamma_fs[k]	= gamma[k] * Fs(delta[k]/t);
		ilt		   += 0.4 / t * creal(gamma_fs[k]);
	}
	
	// Return inverse Laplace transform
	return ilt;
}


