
/***************************************************************************************
*	Inverse Laplace transform using Talbot's algorithm, compared to Euler Inversion    *
*	Inspired by MATLAB code by Tucker McClure                                          *
*	Method: Talbot's method, according to Abate and Whitt (2006)                       *
*                                                                                      *
*	Fredrik Gjertsen 	<fredrik.gjertsen@gmail.com>  		                           *
***************************************************************************************/

#include <stdlib.h>		/* Standard Library functions			*/
#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standard Library of Complex Numbers 	*/
#include <math.h>		/* Standard Library of Math				*/

#define M		64
#define cot(x)	(cos(x)/sin(x))

// Set what Laplace transform function to invert (see list of functions below)
static int nLaplaceTransform;

// Set time vector
static const double ti[] = {0.05, 1.05, 2.05, 3.05, 4.05};


/***************************************************************************************
* Laplace transform function 														   *
***************************************************************************************/
double complex F_s(double complex s)
{	
	switch (nLaplaceTransform) {
	// Laplace transforms to invert:	
	/* Alt 0: Ramp function, type 1		*/	case 0:		return ( 2.0/(s*s) );	
	/* Alt 1: Step function				*/	case 1:		return ( 1.0/s );
	/* Alt 2: Ramp function, type 2		*/	case 2:		return ( 1.0/(s*s) );
	/* Alt 3: Damped sine				*/	default:	return ( 2.0*M_PI/((s+2.0)*(s+2.0) + 4.0*M_PI*M_PI) );
	}
}

/***************************************************************************************
* Binomial coeffient ("from n choose k")                                               *
***************************************************************************************/
int bnml(int n, int k)
{
    int ans = 1;
    k = k > n-k ? n-k : k;
    
    for (int j = 1; j <= k; j++, n--){
        if		(n%j==0)        ans	*= n/j;
        else if (ans%j==0)		ans	 = ans/j*n;
        else    				ans	 = (ans*n)/j;
    }
    return ans;
}

/***************************************************************************************
* Talbot inversion formula                                                             *
***************************************************************************************/
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
	
	// Calculate inverse LT
	double ilt = 0.0;
	for (int k = 0; k < M; k++){
		gamma_fs[k]	= gamma[k] * Fs(delta[k]/t);
		ilt		   += 0.4 / t * creal(gamma_fs[k]);
	}
	
	// Return inverse Laplace transform
	return ilt;
}

/***************************************************************************************
* Euler inversion formula                                                              *
***************************************************************************************/
double Euler_inversion(						// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	)
{
	int m = M/2;
	double complex
		beta[M + 1];
	double
		xi[M + 1],
		eta[M + 1];

	// Calculate xi
	xi[0] = 0.5;
	xi[m] = 1.0;
	xi[M] = 1.0 / pow(2.0, m);
	for (int k = 1; k < m; k++){
		xi[k]		= 1.0;
		xi[M - k]	= xi[M - k + 1] + bnml(m, k) / pow(2.0, m);
	}

	// Calculate beta and eta
	for (int k = 0; k <= M; k++){
		beta[k] = m*log(10.0)/3.0 + M_PI*k*I;
		eta[k]	= (1 - (k%2)*2) * xi[k];
	}
	
	// Calculate inverse LT
	double ilt = 0.0;
	for (int k = 0; k <= M; k++)
		ilt	 += pow(10.0, m/3.0) / t * eta[k] * creal( Fs(beta[k]/t) );
	
	// Return inverse Laplace transform
	return ilt;
}


/***************************************************************************************
* Main entry point                                                                     *
***************************************************************************************/
int main(int argc, char** argv)
{	
	// Allocate resources
	int ni = (int)(sizeof(ti)/sizeof(ti[0]));
	double *ft = (double*)malloc(2 * ni * sizeof(double));
	
	if (argc > 1) {
		nLaplaceTransform = atoi(argv[1]);
	} else {
		nLaplaceTransform = 3;
	}
	
	// Calculate inverse Laplace transforms
	for (int i = 0; i < ni; i++) {
		ft[i] 		= Talbot_inversion( F_s, ti[i] );		// Talbot inversion
		ft[i + ni]	= Euler_inversion(  F_s, ti[i] );		// Euler inversion
	}
	
	// Print results to console
	printf("\nTime\t\tTalbot inversion\tEuler inversion\n");
	for (int i = 0; i < ni; i++) {
		printf("t = %.2f\tf(t) = %.3E\tf(t) = %.3E\n", ti[i], ft[i], ft[i+ni]);
	}
	
	// Free allocated memory
	free(ft);

    return 0;
}


