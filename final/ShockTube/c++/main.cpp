#include <iostream>
#include <cmath>
#include <fstream>

#include "hFiles/analytical.h"
#include "hFiles/Upwind.h"

using namespace std;

void analytical( const int , const double , const double , const double , const double ,
                const double , const double , const double , const double ,
                const double , const double , double * ) ;

void Upwind ( const int , const double , double , int , double , const double ,
             double * , double * , double * , const double , double * ) ;

void Timestep( const double , const double , const double , const double , const double ,
              const double , double * , const double ) ;


int main()
{
    const int nx = 2001 ;        // nodes and set by user
    const double tfinal = 0.849694 ;  // set by user

    //boundaries (can be set)
    const double x_min = 0.0 ;
    const double x_max = 10.0 ;
    const double h = 1.00 * ( x_max - x_min ) / ( nx - 1 ) ;  // ( x_max - x_min ) Tube length ( meters )

    const double gamma = 1.4 ;

    const double cfl = 0.5 ;   // set by user

    const int x0 = 5 ; //diaphram start location ( set by user )
    const int n0 = x0 / h + 1 ;  //node of diaphram star


    //Initial conditions ( set by user )
    const double rho_left = 1.0 ;
    const double p_left = 1.0 ;
    const double u_left = 0.0 ;

    const double rho_right = 0.125 ;
    const double p_right = 0.1 ;
    const double u_right = 0.0 ;


    double time = 0 ;


    double *x , *rho , *rhou , *E ;
    x = new double [nx] ;
    rho = new double [nx] ;
    rhou = new double [nx] ;
    E = new double [nx] ;

    double **Q ;
    Q = new double *[3] ;

    for ( int i = 0 ; i < 3 ; i++ )
        Q[i]    = new double [nx] ;


    int i , maxstep;
    double dt ;


    x[0] = x_min ;
    for ( int i = 1 ; i < nx ; i++ )
        x[i] = x[i-1] + h ;           //Grid

    for ( i = 0 ; i < nx ; i++ )
    {
        rhou[i] = 0.0 ;
        if ( i < n0 )
        {
            rho[i] = rho_left ;
            E[i] = p_left / ( gamma - 1 ) ;
        }
        else
        {
            rho[i] = rho_right ;
            E[i] = p_right / ( gamma - 1 ) ;
        }
        Q[0][i] = rho[i] ;
        Q[1][i] = rhou[i] ;
        Q[2][i] = E[i] ;
    }


    Timestep( gamma , p_right , p_left , rho_right , rho_left , h , &dt , cfl ) ;

    maxstep = tfinal / dt ;
    //cout << maxstep << endl ;
    maxstep += 544 ;

    Upwind ( nx , h , dt , maxstep , time , gamma , rho , rhou , E , cfl , x ) ;

    analytical( nx , tfinal , h , x0 , gamma , p_right , p_left , rho_right , rho_left , u_right , u_left , x ) ;


    delete []x ;
    delete []rho ;
    delete []rhou ;
    delete []E ;

    for ( i = 0 ; i < 3 ; i++ )
        delete []Q[i] ;

    delete []Q ;

    return 0;
}


//============================================================


void Timestep( const double gamma , const double p_right , const double p_left , const double rho_right , const double rho_left ,
              const double h , double *dt , const double cfl )
{
    double duml , dumr , cmax ;
    dumr = gamma * p_right / rho_right ;
    duml = gamma * p_left / rho_left ;

    cmax = sqrt ( duml ) ;
    if ( dumr > duml )
        cmax = sqrt ( dumr ) ;

    *dt = cfl * h / cmax ;
}
