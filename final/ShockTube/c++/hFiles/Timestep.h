#ifndef TIMESTEP_H_INCLUDED
#define TIMESTEP_H_INCLUDED


//Time step
//============================================================


double Timestep( double *c , double *u , const double h , const int nx , const double cfl )
{
    double cmax = 1.0 ;

       for ( int i = 0 ; i < nx ; i++ )
    {
        if ( fabs ( u[i] ) + fabs ( c[i] ) >  cmax )
            cmax = fabs ( u[i] ) + fabs ( c[i] ) ;
    }

    return ( cfl * h / fabs ( cmax ) ) ;
}


//============================================================


#endif // TIMESTEP_H_INCLUDED
