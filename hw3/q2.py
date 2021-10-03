import numpy as np
import matplotlib.pyplot as plt
dxz = [0.2,0.05,0.025,0.006,0.003]

start = 0
finish = 15

name = 'Euler'
for dx in dxz:
    n = int((finish-start)/dx)+1
    x = np.linspace(0,15,n)
    y = np.full_like(x, 1)
    yp = np.full_like(x, 0)
    for i in range(n-1):
        yp[i] = -0.2*y[i]-2*np.cos(2*x[i])*y[i]**2
        y[i+1] = y[i] + dx*yp[i]
    plt.plot(x,y)
y = 1 / (100/101*np.sin(2*x)-10/101*np.cos(2*x)+111/101*np.exp(0.2*x))
plt.plot(x,y)

plt.legend(['0.2','0.05','0.025','0.006','0.003','exact'])
plt.gca().set_ylim([0.0,1.4])
plt.savefig('plot/1-a-'+name+'.png')
plt.title(name)
plt.show()
