#ifndef G_H_INCLUDED
#define G_H_INCLUDED


//============================================================


double G ( double p , const double gamma , const double alpha , const double p_right , const double rho_right )
{
	return ( ( ( 1 - ( 1 / alpha ) ) / sqrt ( rho_right * ( p + p_right / alpha ) ) ) *
         ( 1 - ( p - p_right ) / ( 2 * pow ( ( p + p_right / alpha ) , 1.5 ) ) ) +
             sqrt ( ( gamma - 1 ) / gamma ) * pow ( p , - ( gamma + 1 ) / ( 2 * gamma ) ) ) ;
}


//============================================================


#endif // G_H_INCLUDED
