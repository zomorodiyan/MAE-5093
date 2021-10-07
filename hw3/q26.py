import numpy as np
import matplotlib.pyplot as plt

def a(x):
    return -(x+3)/(x+1)
def b(x):
    return (x+3)/(x+1)**2
def f(x):
    return 2*(x+1) + 3*b(x)
def step(x,s,t,h):
    k1t = h*s
    k1s = h*(-a(x)*s-b(x)*t+f(x))
    k2t = h*(s+k1s/2)
    k2s = h*(-a(x+h/2)*(s+k1s/2)-b(x+h/2)*(t+k1t/2)+f(x+h/2))
    k3t = h*(s+k2s/2)
    k3s = h*(-a(x+h/2)*(s+k2s/2)-b(x+h/2)*(t+k2t/2)+f(x+h/2))
    k4t = h*(s+k3s)
    k4s = h*(-a(x+h)*(s+k3s)-b(x+h)*(t+k3t)+f(x+h))
    return [t + 1/6*k1t + 1/3*(k2t+k3t) + 1/6*k4t\
           ,s + 1/6*k1s + 1/3*(k2s+k3s) + 1/6*k4s]

tA = 5; tB = 4
start = 0; finish = 2; dx = 0.02
n = int((finish-start)/dx)
x = np.linspace(start,finish,n);
t1 = np.empty(n); t2 = np.empty(n); s1 = np.empty(n); s2 = np.empty(n)

t1[0] = tA; t2[0] = tA
s1[0] = 50; s2[0] = -28 #guess
#while(np.abs(s1[0]-s2[0])>1):
#for k in range(5):
while(np.abs(t2[-1]-tB)>1e-6):
    for i in range(n-1):
        t1[i+1],s1[i+1] = step(x[i],s1[i],t1[i],dx)
        t2[i+1],s2[i+1] = step(x[i],s2[i],t2[i],dx)
    guess = s2[0] - (t2[-1]-tB) * (s2[0]-s1[0]) / (t2[-1]-t1[-1])
    s1[0] = s2[0]
    s2[0] = guess

plt.plot(x,t2)
plt.show()
