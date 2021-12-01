from tools import *
import numpy as np

def main():
    nx = 2001 # nodes and set by user
    tfinal = 0.849694  # set by user

    #boundaries (can be set)
    x_min = 0.0
    x_max = 10.0
    h = 1.00 * ( x_max - x_min ) / ( nx - 1 ) # ( x_max - x_min ) Tube length ( meters )
    gamma = 1.4
    cfl = 0.5 # set by user
    x0 = 5 # diaphram start location ( set by user )
    n0 = x0 / h + 1 # node of diaphram star

    #Initial conditions ( set by user )
    rho_left = 1.0
    p_left = 1.0
    u_left = 0.0

    rho_right = 0.125
    p_right = 0.1
    u_right = 0.0

    time = 0

    x = np.empty(nx)
    rho = np.empty(nx)
    rhou = np.empty(nx)
    E = np.empty(nx)
    Q = np.empty((3,nx))

    x = np.linspace(x_min,x_max,nx) #qqq not completely sure if it is correct

    for i in range(nx):
        rhou[i] = 0.0 ;
        if ( i < n0 ):
            rho[i] = rho_left
            E[i] = p_left / ( gamma - 1 )
        else:
            rho[i] = rho_right ;
            E[i] = p_right / ( gamma - 1 )
        Q[0][i] = rho[i]
        Q[1][i] = rhou[i]
        Q[2][i] = E[i]

    dt = timestep2(gamma, p_right, p_left, rho_right, rho_left, h, cfl) ;
    maxstep = tfinal / dt ;
    maxstep += 544
    Upwind( nx , h , dt , maxstep , time , gamma , rho , rhou , E , cfl , x )
    analytical( nx , tfinal , h , x0 , gamma , p_right , p_left , rho_right ,\
            rho_left , u_right , u_left , x ) ;

