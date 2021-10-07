from sympy import symbols, diff, poly, Function, factorial
import numpy as np
x, h, a, b = symbols('x, h, a, b')
f = Function('f')
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))
eq = h**2*(t(-1,5) + b*f(x) + t(1,5)).diff(x,2) -a*(t(1,7)-2*f(x)+t(-1,7))

mat = np.zeros([2,3])
count = 0
for order in range(0,6):
    p = poly(eq,[f(x).diff(x,order)])
    c = p.coeffs()
    if len(c)==2:
        p = poly(c[0],[h])
        c = p.coeffs()[0]
        if(len(poly(c,[a]).coeffs())==2):
            mat[count,0] = poly(c,[a]).coeffs()[0]
        if(len(poly(c,[b]).coeffs())==2):
            mat[count,1] = poly(c,[b]).coeffs()[0]
        mat[count,2] = -1*(c - a*mat[count,0] - b*mat[count,1])
        count += 1
        print(c)
    else:
        print('none')
        pass
l = np.array(mat[:,:2])
r = np.array(mat[:,2])
r = np.expand_dims(np.array(mat[:,2]),axis=1)
print('ans: \n',np.linalg.inv(l)@r)
