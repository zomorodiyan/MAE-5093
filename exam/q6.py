from sympy import symbols,Function,factorial,diff,poly ,zeros,Matrix,shape
import matplotlib.pyplot as plt
import numpy as np

x, h, b2, a1, a2, a3 = symbols('x, h, b2, a1, a2, a3'); f = Function('f')
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))

#f(x) ≡ yb
y1 = t(0.5,5); y1pp = y1.diff(x,2)
y2 = t(1.5,5); y2pp = y2.diff(x,2)
y3p = t(2.5,5)
E = y1pp + b2*y2pp - a1*y1 - a2*y2 - a3*f(x).diff(x,1)

mat1 = zeros(4,4); mat2 = zeros(4,1); vars = Matrix([a1,a2,a3,b2])
for i in range(shape(mat1)[0]):
    p = poly(E,[f(x).diff(x,i)]); c = p.coeffs()[0]
    mat2[i] = -c
    for j in range(shape(mat2)[0]):
        C = poly(c,[vars[j]]).coeffs()
        if len(C)==2:
            mat1[i,j] = C[0]
            mat2[i,0] += C[0]*vars[j]
print(mat1.inv() * mat2)
