#ifndef NEWTONRAPHSON_H_INCLUDED
#define NEWTONRAPHSON_H_INCLUDED

#include "F.h"
#include "G.h"

double F ( double , const double , const double , const double , const double ) ;
double G ( double , const double , const double , const double , const double ) ;


//Newton Raphson solver
//============================================================


void NewtonRaphson ( double p0 , double *x , const double gamma , const double alpha , const double p_right , const double rho_right )
{
	double x1 = 0 , error = 0.001 ;

	*x = p0 - ( F( p0 , gamma , alpha , p_right , rho_right ) ) / ( G( p0 , gamma , alpha , p_right , rho_right ) ) ;

	while ( fabs ( *x - p0 ) > error )
	{
		p0 = *x ;
		*x = p0 - ( F( p0 , gamma , alpha , p_right , rho_right ) ) / ( G( p0 , gamma , alpha , p_right , rho_right ) ) ;
	}
}


//============================================================


#endif // NEWTONRAPHSON_H_INCLUDED
