
/***************************************************************************************
* Examples for inverse Laplace transforms											   *
***************************************************************************************/

#include "nilt.h"		// Numerical inverse Laplace transforms

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
	printf("\n");
	
	// Free allocated memory
	free(ft);

    return 0;
}