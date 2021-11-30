#ifndef F_H_INCLUDED
#define F_H_INCLUDED


//============================================================


double F ( double p , const double gamma , const double alpha , const double p_right , const double rho_right )
{
	return ( ( p - p_right ) * ( 1 - ( 1 / alpha ) ) / sqrt ( ( rho_right * ( p + p_right / alpha ) ) )
         - 2 * ( sqrt ( gamma ) / ( gamma - 1 ) ) * ( 1 - pow( p , ( gamma - 1 ) / ( 2 * gamma ) ) ) ) ;
}


//============================================================


#endif // F_H_INCLUDED
