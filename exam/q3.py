from sympy import symbols, Function, factorial, diff, poly
import matplotlib.pyplot as plt
import numpy as np
x, h, a, b = symbols('x, h, a, b')
f = Function('f')
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))
print("---------part-a----------")
#define equation
eq = (t(1,3) - f(x))/h
#test
p = poly(eq,[f(x).diff(x,1)]); c = p.coeffs()[0]
print("coeff for dy/dt should be 1 and is:",c)
#answer
p = poly(eq,[f(x).diff(x,2)]); c = p.coeffs()[0]
print("the truncation error is",str(c)+'*d2y/dt2')
print("so the order of error is one")
print("---------part-b----------")
delta = 0.025
x = np.arange(-2.5, 0.5, delta)
y = np.arange(-1.5, 1.5, delta)
X, Y = np.meshgrid(x, y)
Z = np.abs(1 + X + Y*1j)
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, [1.0], colors=('black'))
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('Sigma')
ax.xaxis.grid(True, zorder=0)
ax.yaxis.grid(True, zorder=0)
ax.set_xlabel('h*lambda_real')
ax.set_ylabel('h*lambda_imag')
plt.show()
print("showing contour for simga")
print("---------part-c----------")
print("stability region concludes: -2 < h*lambda_real < 0,")
print("if lamda_real==-20 then: 0 < h < 0.1")

