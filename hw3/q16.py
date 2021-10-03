import numpy as np
import matplotlib.pyplot as plt

sigma = 10.0; b = 8/3 # r: control parameter, r>24.74 => kaotic
#ODEz: x'=sigma*(y-x), y'=rx-y-xz, z'=xy-bz

def plot(par1,par2,par3,name1,name2,name3,part):
    plt.plot(x,y); plt.plot(x,z)
    plt.legend(['y', 'z'])
    plt.title(name1+' and '+name2+' for differnet '+name3)
    plt.savefig('plot/'+part+'-'+name1+name2+'per'+name3)
    #plt.show()
    print(name1+name2+'.png is saved in ./plot')
    plt.clf()
def plotc(par1,par2,t,name):
    plt.plot(t,par1); plt.plot(t,par2); plt.legend([name+'1',name+'2'])
    plt.title(name+' evolutions in time')
    plt.savefig('plot/16-c-'+name)
    plt.show()
    print('16-c'+name+'.png is saved in ./plot')
    plt.clf()

def step(xi,yi,zi,dt):
    k1x = dt*sigma*(yi-xi)
    k1y = dt*(r*xi-yi-xi*zi)
    k1z = dt*(xi*yi-b*zi)
    k2x = dt*sigma*((yi+k1y/2)-(xi+k1x/2))
    k2y = dt*(r*(xi+k1x/2)-(yi+k1y/2)-(xi+k1x/2)*(zi+k1z/2))
    k2z = dt*((xi+k1x/2)*(yi+k1y/2)-b*(zi+k1z/2))
    k3x = dt*sigma*((yi+k2y/2)-(xi+k2x/2))
    k3y = dt*(r*(xi+k2x/2)-(yi+k2y/2)-(xi+k2x/2)*(zi+k2z/2))
    k3z = dt*((xi+k2x/2)*(yi+k2y/2)-b*(zi+k2z/2))
    k4x = dt*sigma*((yi+k3y)-(xi+k3x))
    k4y = dt*(r*(xi+k3x)-(yi+k3y)-(xi+k3x)*(zi+k3z))
    k4z = dt*((xi+k3x)*(yi+k3y)-b*(zi+k3z))
    return [xi + 1/6*k1x + 1/3*(k2x+k3x) + 1/6*k4x\
           ,yi + 1/6*k1y + 1/3*(k2y+k3y) + 1/6*k4y\
           ,zi + 1/6*k1z + 1/3*(k2z+k3z) + 1/6*k4z]

start = 0; finish = 25
dt = 0.1
n = int((finish-start)/dt)+1
x = np.empty(n); y = np.empty(n); z = np.empty(n)
t = np.linspace(0,25,n)
#(a) 0<=t<=25
r = 20; x[0],y[0],z[0]= [1,1,1]

for i in range(n-1):
    x[i+1],y[i+1],z[i+1] = step(x[i],y[i],z[i],dt)
plot(x,y,z,'x','y','z','16-a')
plot(x,z,y,'x','z','y','16-a')
plot(y,z,x,'y','z','x','16-a')

#(b)
r = 27.00; x[0],y[0],z[0]= [0.01,0.01,0.01]
for i in range(n-1):
    x[i+1],y[i+1],z[i+1] = step(x[i],y[i],z[i],dt)
plot(x,y,z,'x','y','z','16-b')
plot(x,z,y,'x','z','y','16-b')
plot(y,z,x,'y','z','x','16-b')

#(c) compare trajectory in time for two nearby initial conditions
r=28
x[0],y[0],z[0] = [6,6,6]
for i in range(n-1):
    x[i+1],y[i+1],z[i+1] = step(x[i],y[i],z[i],dt)
x2 = np.copy(x); y2 = np.copy(y); z2 = np.copy(z)

x2[0],y2[0],z2[0] = [6,6.01,6]
for i in range(n-1):
    x2[i+1],y2[i+1],z2[i+1] = step(x2[i],y2[i],z2[i],dt)
plotc(x,x2,t,'x'); plotc(y,y2,t,'y'); plotc(z,z2,t,'z')



