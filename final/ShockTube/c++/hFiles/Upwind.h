#ifndef UPWIND_H_INCLUDED
#define UPWIND_H_INCLUDED
#include "Timestep.h"

//============================================================


void Upwind ( const int nx , const double h , double dt , int maxstep , double time , const double gamma ,
             double *rho , double *rhou , double *E , const double cfl , double *x )
{
    ofstream fDensity ("UpwindDensity.plt") ;
    fDensity.setf     ( ios :: showpoint ) ;

    ofstream fVelocity ("UpwindVelocity.plt") ;
    fVelocity.setf     ( ios :: showpoint ) ;

    ofstream fPressure ("UpwindPressure.plt") ;
    fPressure.setf     ( ios :: showpoint ) ;

    ofstream fInternalEnergy ("UpwindInternalEnergy.plt") ;
    fInternalEnergy.setf     ( ios :: showpoint ) ;

    ofstream fMakh ("UpwindMakh.plt") ;
    fMakh.setf     ( ios :: showpoint ) ;

    ofstream fMssFlow ("UpwindMassFlow.plt") ;
    fMssFlow.setf     ( ios :: showpoint ) ;

    ofstream fEntropy ("UpwindEntropy.plt") ;
    fEntropy.setf     ( ios :: showpoint ) ;

    double *p , *c , *u , *m , *F1 , *F2 , *F3 ;
    p = new double [nx] ;
    c = new double [nx] ;
    u = new double [nx] ;
    m = new double [nx] ;
    F1 = new double [nx] ;
    F2 = new double [nx] ;
    F3 = new double [nx] ;

    const double Cv = 0.7179 ;

    for ( int istep = 0 ; istep < maxstep ; istep++ )
    {
        if ( istep == maxstep - 1 )
        {
            fDensity << "Variables=X,Denstity\nZONE T=\"Numerical" << "\"\nSOLUTIONTIME= " << time << endl ;
            fVelocity << "Variables=X,Velocity\nZONE T=\"Numerical" << "\"\nSOLUTIONTIME= " << time << endl ;
            fPressure << "Variables=X,Pressure\nZONE T=\"Numerical" << "\"\nSOLUTIONTIME= " << time << endl ;
            fInternalEnergy << "Variables=X,Internal Energy\nZONE T=\"Numerical" << time << "\"\nSOLUTIONTime= " << time << endl ;
            fMakh << "Variables=X,Mach\nZONE T=\"Numerical" << "\"\nSOLUTIONTime= " << time << endl ;
            fMssFlow << "Variables=X,Mass Flow\nZONE T=\"Numerical" << "\"\nSOLUTIONTime= " << time << endl ;
            fEntropy << "Variables=X,Entropy\nZONE T=\"Numerical" << "\"\nSOLUTIONTime= " << time << endl ;
        }



        for ( int i = 0 ; i < nx ; i++ )
        {
            p[i] = ( gamma - 1 ) * ( E[i] - 0.5 * ( rhou[i] * rhou[i] / rho[i] ) ) ;
            c[i] = sqrt( gamma * p[i] / rho[i] ) ;
            u[i] = rhou[i] / rho[i] ;
            m[i] = u[i] / c[i] ;
        }

        dt = Timestep( c , u , h , nx , cfl ) ;

        for ( int i = 1 ; i < nx - 1 ; i++ ) // Find fluxes
        {
            F1[i] = 0.5 * ( rhou[i+1] + rhou[i] ) - 0.5 * ( fabs( rhou[i+1] ) - fabs( rhou[i] ) ) ;

            F2[i]= 0.5 * ( u[i+1] * rhou[i+1] + p[i+1] + u[i] * rhou[i] + p[i] )
            - 0.5 * ( fabs( u[i+1] ) * rhou[i+1] - fabs( u[i] ) * rhou[i] )
            - 0.5 * ( p[i+1] * m[i+1] - p[i] * m[i] ) ;

            F3[i] = 0.5 * ( u[i+1] * ( E[i+1] + p[i+1] ) + u[i] * ( E[i] + p[i] ) )
            - 0.5 * ( fabs( u[i+1] ) * E[i+1] - fabs( u[i] ) * E[i] )
            - 0.5 * ( p[i+1] * c[i+1] - p[i] * c[i] ) ;

            if ( m[i] > 1)
            {
                F2[i] = rhou[i] * u[i] +p[i] ;
                F3[i] = ( E[i] + p[i] ) *u[i] ;
            }


            if ( m[i] < -1 )
            {
                F2[i] = rhou[i+1] * u[i+1] + p[i+1] ;
                F3[i] = ( E[i+1] + p[i+1] ) * u[i+1] ;
            }
        }

        for ( int i = 1 ; i < nx - 1 ; i++ ) // Update solution
        {
            rho[i] = rho[i] - ( dt / h ) * ( F1[i] - F1[i-1] ) ;
            rhou[i] = rhou[i] - ( dt / h ) * ( F2[i] - F2[i-1] ) ;
            E[i] = E[i] - ( dt / h ) * ( F3[i] - F3[i-1] ) ;
        }

        for ( int i = 0 ; i < nx ; i++)
        {
            if ( istep == maxstep - 1 )
            {
                p[i] = ( ( gamma - 1 ) * ( E[i] - 0.5 * rho[i] * u[i] * u[i] ) ) ;
                m[i] = u[i] / sqrt ( gamma * p[i] / rho[i] ) ;
                fDensity << x[i] << '\t' << rho[i] << endl ;
                fVelocity << x[i] << '\t' << rhou[i] / rho[i] << endl ;
                fPressure << x[i] << '\t' << p[i] << endl ;
                fInternalEnergy << x[i] << '\t' << E[i] << endl ;
                fMakh << x[i] << '\t' << rhou[i] / ( sqrt ( gamma * p[i] / rho[i] ) * rho[i] ) << endl ;
                fMssFlow << x[i] << '\t' << rho[i] * u[i] << endl ;
                fEntropy << x[i] << '\t' << Cv * log ( p[i] / pow ( rho[i] , gamma ) ) << endl ;
            }

        }

        time += dt ;

    }

    fDensity.close() ;
    fVelocity.close() ;
    fPressure.close() ;
    fInternalEnergy.close() ;
    fMakh.close() ;
    fMssFlow.close() ;
    fEntropy.close() ;

    delete []p ;
    delete []c ;
    delete []u ;
    delete []m ;
    delete []F1 ;
    delete []F2 ;
    delete []F3 ;
}


//============================================================


#endif // UPWIND_H_INCLUDED
