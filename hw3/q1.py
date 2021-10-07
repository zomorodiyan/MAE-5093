import numpy as np
import matplotlib.pyplot as plt

#ODE: y'=-(2+0.01*x**2)*y
#exact solution: y = (1 + 3*x + 3*x**2 + x**3) / np.exp(x)

def alpha(xi):
    return 3*xi/(1+xi)
def beta(xi):
    return 2*(1+xi)**3*np.exp(-xi)
def step(name,xi,yi,dx):
    if name=='Euler':
        yp = -(2+0.01*xi**2)*yi
        return yi + dx*yp
    if name=='BackEuler':
        return yi/(1+dx*(2+0.01*xi**2))
    if name=='Trapezoidal':
        return (yi - dx/2*(2+0.01*xi**2)*yi) / (1+dx/2*(2+0.01*(xi+dx)**2))
    if name=='RK2':
        yhalf = yi - dx/2*(2+0.01*xi**2)*yi
        return yi - dx*(2+0.01*(xi+dx/2)**2)*yhalf
    if name=='RK4':
        k1 = -dx*(2+0.01*xi**2)*yi
        k2 = -dx*(2+0.01*(xi+dx/2)**2)*(yi+k1/2)
        k3 = -dx*(2+0.01*(xi+dx/2)**2)*(yi+k2/2)
        k4 = -dx*(2+0.01*(xi+dx)**2)*(yi+k3)
        return y[i] + 1/6*k1 + 1/3*(k2+k3) + 1/6*k4

dxz = [1.0,0.5,0.1]
names = ['Euler','BackEuler','Trapezoidal','RK2','RK4']
start = 0; finish = 15

for name in names:
    for dx in dxz:
        n = int((finish-start)/dx)+1
        x = np.linspace(0,15,n)
        y = np.full_like(x, 1)
        yp = np.full_like(x, 0)
        for i in range(n-1):
            y[i+1] = step(name,x[i],y[i],dx)
        plt.plot(x,y)
    plt.legend(['dx=1.0', 'dx=0.5', 'dx=0.1'])
    plt.gca().set_ylim([-0.25,1.25])
    plt.title(name)
    plt.savefig('plot/1-'+name+'.png')
    #plt.show()
    plt.clf()
    print('1-'+name+'.png is saved in ./plot')
