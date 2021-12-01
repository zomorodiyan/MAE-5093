import numpy as np

def F(pressure,gamma,alpha,pressure_r,density_r):
    return (pressure-pressure_r)*(1-(1/alpha))\
        /np.sqrt((density_r*(pressure+pressure_r/alpha)))\
        -2*(np.sqrt(gamma)/(gamma-1))*(1-pow(pressure,(gamma-1)/(2*gamma)))

def G(pressure, gamma, alpha, pressure_r, density_r):
    return ((1-(1/alpha))/np.sqrt(density_r*(pressure+pressure_r/alpha)))*\
        (1-(pressure-pressure_r)/(2*pow((pressure+pressure_r/alpha),1.5)))+\
        np.sqrt((gamma-1)/gamma)*pow(pressure,-(gamma+1)/(2*gamma))

def NewtonRaphson(pressure_0, gamma, alpha, pressure_r, density_r):
    # updates x (originally by reference but now returns new value)
    error=0.001
    x=pressure_0 - F(pressure_0, gamma, alpha, pressure_r, density_r)/\
            G(pressure_0, gamma, alpha, pressure_r, density_r)
    while(np.abs(x-pressure_0) > error):
        pressure_0 = x
        x = pressure_0 - F(pressure_0, gamma, alpha, pressure_r, density_r)/\
            G(pressure_0, gamma, alpha, pressure_r, density_r)
    return x

def timestep(c,u,h,nx,cfl):
    cmax = 1.0
    for i in range(nx):
        if (np.abs(u[i]) + np.abs(c[i]) > cmax):
            cmax = np.abs(u[i]) + np.abs(c[i])
    return cfl*h/np.abs(cmax)


def timestep2(gamma, p_right, p_left, rho_right,rho_left, h, cfl):
    dumr = gamma * p_right / rho_right
    duml = gamma * p_left / rho_left
    cmax = np.sqrt (duml) ;
    if dumr > duml:
        cmax = np.sqrt(dumr)
    return cfl * h / cmax

def analytical(nx,t,h,x0,gamma,pressure_r,pressure_l,density_r,density_l,\
        velocity_r,velocity_l,x):
    alpha = (gamma+1)/(gamma-1)

    #speed of sound
    c_l = np.sqrt(gamma*pressure_l/density_l)
    c_r = np.sqrt(gamma*pressure_r/density_r)

    density = np.empty((nx))
    pressure = np.empty((nx))
    velocity = np.empty((nx))
    energy = np.empty((nx))

    pressure_0 = 0.1 ;
    pressure_post = NewtonRaphson(pressure_0, gamma, alpha, pressure_r, density_r)
    velocity_post = 2*(np.sqrt(gamma)/(gamma-1))*(1-pressure_post**((gamma-1)/(2*gamma)))
    density_post = density_r*((pressure_post/pressure_r)+(1/alpha))/\
            (1+(pressure_post/pressure_r)/alpha)
    velocity_shock = velocity_post*((density_post/density_r)/((density_post/density_r)-1))
    density_middle = density_l*(pressure_post/pressure_l )**1/gamma

    #Key Positions
    x1 = x0 - c_l * t ;
    x3 = x0 + velocity_post * t ;
    x4 = x0 + velocity_shock * t ;

    #determining x2
    c_2 = c_l-(gamma-1)/2*velocity_post;
    x2 = x0+(velocity_post-c_2)*t

    Cv = 0.7179 ;

    for i in range(nx):
        if x[i]<x1:
            #Solution before x1
            density[i] = density_l
            pressure[i] = pressure_l
            velocity[i] = velocity_l
        elif x[i] < x2:
            #Solution between x1 and x2
            c = (x0-x[i])/t/alpha+(1-1/alpha)*c_l
            density[i] = density_l*c/c_l**2/(gamma-1)
            pressure[i] = pressure_l*(density[i]/density_l)**gamma
            velocity[i] = (1-1/alpha)*(-1*(x0-x[i])/t + c_l)
        elif x[i] < x3:
            #Solution between x2 and x3
            density[i] = density_middle
            pressure[i] = pressure_post
            velocity[i] = velocity_post
        elif x[i] < x4:
            #Solution b/w x3 and x4
            density[i] = density_post ;
            pressure[i] = pressure_post ;
            velocity[i] = velocity_post ;
        elif x[i] > x4:
            #Solution after x4
            density[i] = density_r ;
            pressure[i] = pressure_r ;
            velocity[i] = velocity_r ;
        else: raise ValueError("x[i] is invalid")

        energy[i] = pressure[i]/(gamma-1)*density[i]

    return density,pressure,velocity,energy



def Upwind(nx, h, dt, maxstep, time, gamma, rho, rhou, E, cfl, x):
    p = np.empty((nx))
    c = np.empty((nx))
    u = np.empty((nx))
    m = np.empty((nx))
    M = np.empty((nx)) # Mach
    Ent = np.empty((nx)) # Entropy
    F1 = np.empty((nx))
    F2 = np.empty((nx))
    F3 = np.empty((nx))

    Cv = 0.7179 ;
    for istep in range(maxstep):
        for i in range(nx):
            p[i] = (gamma-1)*(E[i]-0.5*(rhou[i]*rhou[i]/rho[i]))
            c[i] = np.sqrt(gamma*p[i]/rho[i])
            u[i] = rhou[i]/rho[i]
            m[i] = u[i]/c[i]
        dt = timestep(c, u, h, nx, cfl)

        for i in range(1,nx-1): # Find fluxes
            F1[i] = 0.5*(rhou[i+1]+rhou[i])-0.5*(np.abs(rhou[i+1])-np.abs(rhou[i]))
            F2[i] = 0.5*(u[i+1]*rhou[i+1]+p[i+1]+u[i]*rhou[i]+p[i])\
            - 0.5*(np.abs(u[i+1])*rhou[i+1]-np.abs(u[i])*rhou[i])\
            - 0.5*(p[i+1]*m[i+1]-p[i]*m[i])

            F3[i] = 0.5*(u[i+1]*(E[i+1]+p[i+1])+u[i]*(E[i]+p[i]))\
            - 0.5*(np.abs(u[i+1])*E[i+1]-np.abs(u[i])*E[i])
            - 0.5*(p[i+1]*c[i+1]-p[i]*c[i])

            if m[i] > 1:
                F2[i] = rhou[i] * u[i] +p[i]
                F3[i] = ( E[i] + p[i] ) *u[i]
            if m[i] < -1:
                F2[i] = rhou[i+1] * u[i+1] + p[i+1]
                F3[i] = ( E[i+1] + p[i+1] ) * u[i+1]
        for i in range(1,nx-1): # Update solution
            rho[i] = rho[i]-(dt/h)*(F1[i]-F1[i-1])
            rhou[i] = rhou[i]-(dt/h)*(F2[i]-F2[i-1])
            E[i] = E[i]-(dt/h)*(F3[i]-F3[i-1])

        time += dt ;

    for i in range(nx):
        p[i] = ((gamma - 1 ) * ( E[i] - 0.5 * rho[i] * u[i] * u[i] ) )
        m[i] = u[i] / np.sqrt ( gamma * p[i] / rho[i] )
        m[i] = rhou[i] / ( np.sqrt ( gamma * p[i] / rho[i] ) * rho[i] )
        Ent[i] = Cv * np.log ( p[i] / pow ( rho[i] , gamma ) )

    return rho,rhou,E

def main():
    print(F(2,2,2,2,2))
    print(G(2,2,2,2,2))
    print(NewtonRaphson(2,2,2,2,2))
    print(timestep([1,2],[1,2],2,2,2))
    print(timestep2(2, 4, 2, 3,1, 1, 0.1))
    d,p,v,e = analytical(6,1,0.1,0.5,2,3,1,4,2,0,0,[0,0.2,0.4,0.6,0.8,1.0]); print(d)
    r,ru,e = Upwind(6, 0.1, 0.1, 10, 1, 2 , [1,1,1,3,3,3],\
            [1,1,1,3,3,3], [1,1,1,2,2,2], 0.03, [0,0.2,0.4,0.6,0.8,1.0])
    print(r)

if __name__ == "__main__":
    main()
