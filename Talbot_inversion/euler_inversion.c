
/***************************************************************************************
*	Inverse Laplace transform, Eulers Inversion Method								   *
*	Inspired by MATLAB code by Tucker McClure                                          *
*                                                                                      *
*	Fredrik Gjertsen 	<fredrik.gjertsen@gmail.com>  		                           *
***************************************************************************************/

#include "euler_inversion.h"

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
* Special power function where the exponent is an integer                              *
***************************************************************************************/
double powi(double dBase, int nExponent)
{
	double val = 1.0;
	for (int i = 0; i < nExponent; i++){
		val *= dBase;
	}	
	return val;
}

/***************************************************************************************
* Euler inversion formula                                                              *
***************************************************************************************/
double Euler_inversion(						// Ret: Inverse Laplace transform of function Fs at t
	double complex (*Fs)(double complex),	// In:	Pointer to Laplace transform function
	double t								// In:	Time at which to evaluate ILT
	)
{
	double complex
		beta[M + 1];
	double
		xi[M + 1],
		eta[M + 1];

	// Calculate xi
	xi[0] = 0.5;
	xi[m] = 1.0;
	xi[M] = 1.0 / powi(2, m);
	for (int k = 1; k < m; k++){
		xi[k]		= 1.0;
		xi[M - k]	= xi[M - k + 1] + bnml(m, k) / powi(2, m);
	}

	// Calculate beta and eta
	for (int k = 0; k <= M; k++){
		beta[k] = m*log(10.0)/3.0 + M_PI*k*I;
		eta[k]	= (1 - (k%2)*2) * xi[k];
	}
	
	// Calculate inverse Laplace Transform
	double ilt = 0.0;
	for (int k = 0; k <= M; k++)
		ilt	 += pow(10.0, m/3.0) / t * eta[k] * creal( Fs(beta[k]/t) );
	
	// Return inverse Laplace transform
	return ilt;
}