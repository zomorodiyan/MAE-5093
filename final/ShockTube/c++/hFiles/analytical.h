#ifndef ANALYTICAL_H_INCLUDED
#define ANALYTICAL_H_INCLUDED

#include "NewtonRaphson.h"

using namespace std ;

void NewtonRaphson ( double , double * , const double , const double , const double , const double ) ;

//Analytical solution
//============================================================


void analytical( const int nx , const double t , const double h , const double x0 , const double gamma ,
                const double p_right , const double p_left , const double rho_right , const double rho_left ,
                const double u_right , const double u_left , double *x )
{
    ofstream fDensity ("AnalayDensity.plt") ;
    fDensity.setf     ( ios :: showpoint ) ;

    ofstream fVelocity ("AnalayVelocity.plt") ;
    fVelocity.setf     ( ios :: showpoint ) ;

    ofstream fPressure ("AnalayPressure.plt") ;
    fPressure.setf     ( ios :: showpoint ) ;

    ofstream fInternalEnergy ("AnalayInternalEnergy.plt") ;
    fInternalEnergy.setf     ( ios :: showpoint ) ;

    ofstream fMakh ("AnalayMakh.plt") ;
    fMakh.setf     ( ios :: showpoint ) ;

    ofstream fMssFlow ("AnalayMassFlow.plt") ;
    fMssFlow.setf     ( ios :: showpoint ) ;

    ofstream fEntropy ("AnalayEntropy.plt") ;
    fEntropy.setf     ( ios :: showpoint ) ;


    const double alpha = ( gamma + 1 ) / ( gamma - 1 ) ;

    //speed of sound
    const double c_l = sqrt( gamma * p_left  / rho_left  ) ;
    const double c_r = sqrt( gamma * p_right / rho_right ) ;

    double *rho , *p , *u , *e ;
    rho = new double [nx] ;   //density
    p   = new double [nx] ;   //pressure
    u   = new double [nx] ;   //velocity
    e   = new double [nx] ;   //internal energy

    double v_post , v_shock , rho_post , p_post , rho_middle ;

    double p0 = 0.1 ;
    NewtonRaphson ( p0 , & p_post , gamma , alpha , p_right , rho_right ) ;

    v_post = 2 * ( sqrt ( gamma ) / ( gamma - 1 ) ) * ( 1 - pow ( p_post , ( gamma - 1 ) / ( 2 * gamma ) ) ) ;

    rho_post = rho_right * ( ( ( p_post / p_right ) + ( 1 / alpha ) ) / ( 1 + ( p_post / p_right ) / alpha ) ) ;

    v_shock = v_post * ( ( rho_post / rho_right ) / ( ( rho_post / rho_right ) - 1 ) ) ;

    rho_middle = ( rho_left ) * pow ( ( p_post / p_left ) , 1 / gamma ) ;


    //Key Positions
    const double x1 = x0 - c_l * t ;
    const double x3 = x0 + v_post * t ;
    const double x4 = x0 + v_shock * t ;



    //determining x2
    double c_2 = c_l - ( ( gamma - 1 ) / 2 ) * v_post ;
    double x2 = x0 + ( v_post - c_2 ) * t ;

    double c ;
    const double Cv = 0.7179 ;

    fDensity << "Variables=X,Denstity\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fVelocity << "Variables=X,Velocity\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fPressure << "Variables=X,Pressure\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fInternalEnergy << "Variables=X,Internal Energy\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fMakh << "Variables=X,Mach\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fMssFlow << "Variables=X,Mass Flow\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;
    fEntropy << "Variables=X,Entropy\nZONE T=\"Exact" << "\"\nSOLUTIONTime= " << t << endl ;

    for ( int index = 0 ; index < nx ; index++ )
    {
        if ( x[index] < x1 )
        {
            //Solution b4 x1
            rho[index] = rho_left ;
            p[index] = p_left ;
            u[index] = u_left ;
        }

        else if ( x1 <= x[index] && x[index] <= x2 )
        {
            //Solution b/w x1 and x2
            c = ( ( ( x0 - x[index] ) / t ) / alpha ) + ( 1 - ( 1 / alpha ) ) * c_l ;
            rho[index] = rho_left * pow ( ( c / c_l ) , 2 / ( gamma - 1 ) ) ;
            p[index] = p_left * pow ( ( rho[index] / rho_left ) , gamma ) ;
            u[index] = ( 1 - ( 1 / alpha ) ) * ( ( - ( x0 - x[index] ) / t ) + c_l ) ;
        }
        else if ( x2 <= x[index] && x[index] <= x3 )
        {
            //Solution b/w x2 and x3
            rho[index] = rho_middle ;
            p[index] = p_post ;
            u[index] = v_post ;
        }
        else if ( x3 <= x[index] && x[index] <= x4 )
        {
            //Solution b/w x3 and x4
            rho[index] = rho_post ;
            p[index] = p_post ;
            u[index] = v_post ;
        }
        else if ( x4 < x[index] )
        {
            //Solution after x4
            rho[index] = rho_right ;
            p[index] = p_right ;
            u[index] = u_right ;
        }

        e[index] = p[index] / ( ( gamma - 1 ) * rho[index] ) ;

        fDensity << x[index] << '\t' << rho[index] << endl ;
        fVelocity << x[index] << '\t' << u[index] << endl ;
        fPressure << x[index] << '\t' << p[index] << endl ;
        fInternalEnergy << x[index] << '\t' << e[index] << endl ;
        fMakh << x[index] << '\t' << u[index] / sqrt ( gamma * p[index] / rho[index] ) << endl ;
        fMssFlow << x[index] << '\t' << rho[index] * u[index] << endl ;
        fEntropy << x[index] << '\t' << Cv * log ( p[index] / pow ( rho[index] , gamma ) ) << endl ;
        //fEntropy << x[index] << '\t' << u[index] / sqrt ( gamma * p[index] / rho[index] ) << endl ;
        //cout << index << ". analy mach = " << u[index] / sqrt ( gamma * p[index] / rho[index] ) << endl ;
    }

    delete []rho ;
    delete []p ;
    delete []u ;
    delete []e ;

    fDensity.close() ;
    fVelocity.close() ;
    fPressure.close() ;
    fInternalEnergy.close() ;
    fMakh.close() ;
    fMssFlow.close() ;
    fEntropy.close() ;
}


//============================================================


#endif // ANALYTICAL_H_INCLUDED
