import numpy as np

def F(pressure,gamma,alpha,pressure_r,density_right):
    return (pressure-pressure_r)*(1-(1/alpha))\
        /np.sqrt((density_right*(pressure+pressure_r/alpha)))\
        -2*(np.sqrt(gamma)/(gamma-1))*(1-pow(pressure,(gamma-1)/(2*gamma)))

def G(pressure, gamma, alpha, pressure_r, density_r):
    return ((1-(1/alpha))/np.sqrt(density_r*(pressure+pressure_r/alpha)))*\
        (1-(pressure-pressure_r)/(2*pow((pressure+pressure_r/alpha),1.5)))+\
        np.sqrt((gamma-1)/gamma)*pow(pressure,-(gamma+1)/(2*gamma))

print(F(2,2,2,2,2))
print(G(2,2,2,2,2))
