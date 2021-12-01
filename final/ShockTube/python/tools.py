import numpy as np

def F(pressure,gamma,alpha,pressure_r,density_right):
    return (pressure-pressure_r)*(1-(1/alpha))\
        /np.sqrt((density_right*(pressure+pressure_r/alpha)))\
        -2*(np.sqrt(gamma)/(gamma-1))*(1-pow(pressure,(gamma-1)/(2*gamma)))

def G(pressure, gamma, alpha, pressure_r, density_r):
    return ((1-(1/alpha))/np.sqrt(density_r*(pressure+pressure_r/alpha)))*\
        (1-(pressure-pressure_r)/(2*pow((pressure+pressure_r/alpha),1.5)))+\
        np.sqrt((gamma-1)/gamma)*pow(pressure,-(gamma+1)/(2*gamma))

def NewtonRaphson(pressure_0, x, gamma, alpha, pressure_r, density_r):
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

def main():
    print(F(2,2,2,2,2))
    print(G(2,2,2,2,2))
    print(NewtonRaphson(2,2,2,2,2,2))
    print(timestep([1,2],[1,2],2,2,2))
if __name__ == "__main__":
    main()
