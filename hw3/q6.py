import numpy as np
import matplotlib.pyplot as plt

#ODE: yp=-alpha*y+beta | alpha = 3*t/(1+t) , beta = 2*(1+t)^3*np.exp(-t)
#exact solution: y = (1 + 3*x + 3*x**2 + x**3) / np.exp(x)

def alpha(xi):
    return 3*xi/(1+xi)
def beta(xi):
    return 2*(1+xi)**3*np.exp(-xi)
def step(name,xi,yi,dx):
    if name=='Euler':
        yp=-alpha(xi)*yi + beta(xi)
        return yi+dx*yp
    if name=='BackEuler':
        return (yi+beta(xi+dx)*dx)/(1+dx*alpha(xi+dx))
    if name=='Trapezoidal':
        return (yi + dx/2*(-alpha(xi)*yi+beta(xi)+beta(xi+dx)))\
                        /(1+dx/2*alpha(xi+dx))
    if name=='RK2':
        yhalf = yi+dx/2*(-alpha(xi)*yi+beta(xi))
        return y[i] + dx*(-alpha(xi+dx/2)*yhalf+beta(xi+dx/2))
    if name=='RK4':
        k1 = dx*(-alpha(xi)*yi+beta(xi))
        k2 = dx*(-alpha(xi+dx/2)*(yi+k1/2)+beta(xi+dx/2))
        k3 = dx*(-alpha(xi+dx/2)*(yi+k2/2)+beta(xi+dx/2))
        k4 = dx*(-alpha(xi+dx)*(yi+k3)+beta(xi+dx))
        return  yi + 1/6*k1 + 1/3*(k2+k3) + 1/6*k4

dxz = [1.1,0.8,0.2]
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
    y = (1 + 3*x + 3*x**2 + x**3) / np.exp(x)
    plt.plot(x,y)
    plt.legend(['dx=1.1', 'dx=0.8', 'dx=0.2','exact'])
    plt.gca().set_ylim([-0.5,4.5])
    plt.title(name)
    plt.savefig('plot/6-'+name+'.png')
    #plt.show()
    plt.clf()
    print('6-'+name+'.png is saved in ./plot')

