from sympy import symbols, Function, factorial, diff, poly
import matplotlib.pyplot as plt
import numpy as np
x, h, a, b = symbols('x, h, a, b')
f = Function('f')
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))
print("---------part-a----------")
#define equation
eq = (t(1,4) - t(-1,4))/(2*h)
#test
p = poly(eq,[f(x).diff(x,1)]); c = p.coeffs()[0]
print("coeff for dy/dt should be 1 and is:",c)
#answer
p = poly(eq,[f(x).diff(x,3)]); c = p.coeffs()[0]
print("the truncation error is",str(c)+'*d3y/dt3')
print("so the order of error is two")
print("---------part-b----------")
delta = 0.025
x = np.arange(-1, 1, delta)
y = np.arange(-4, 4, delta)
X, Y = np.meshgrid(x, y)
Z = np.abs(X+Y*1j-((X+Y*1j)**2+1)**0.5)
Z2 = np.abs(X+Y*1j+((X+Y*1j)**2+1)**0.5)
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, [0.8,0.9, 1.0], colors=('black','black','black','black'))
CS = ax.contour(X, Y, Z2, [0.8,0.9, 1.0], colors=('black','black','black','black'))
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('Sigma')
ax.xaxis.grid(True, zorder=0)
ax.yaxis.grid(True, zorder=0)
ax.set_xlabel('h*lambda_real')
ax.set_ylabel('h*lambda_imag')
plt.show()
print("showing contour for simga")
print("---------part-c----------")
print("the method is explicit since a direct computation in terms of known \
quntities can be made to get the new values and no iterative technique is \
needed")
