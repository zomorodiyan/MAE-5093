from sympy import symbols,Function,factorial,diff,poly ,zeros,Matrix,shape
import matplotlib.pyplot as plt
import numpy as np

def check_derivation(a1,a2,a3,a4,a5,a6):
    if (a1==-2 and a2==15 and a3==-60 and a4==20 and a5==30 and a6==-3):
        print("the eqation is derived from taylor series definition")
    else:
        print("derivation failed, coeffs don't match to the given ones")
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))

x, h, a1, a2, a3, a4, a5, a6 = symbols('x, h, a1, a2, a3, a4, a5, a6')
f = Function('f')

order=6
E = (a1*t(-3,order)+a2*t(-2,order)+a3*t(-1,order)+a4*f(x)+a5*t(1,order)\
        +a6*t(2,order))/(60*h) -f(x).diff(x,1)

vars = Matrix([a1,a2,a3,a4,a5,a6]); nvars = shape(vars)[0]
mat1 = zeros(nvars,nvars)
mat2 = zeros(nvars,1)
eq_count = 0
for i in range(order):
    if eq_count == nvars: raise ValueError('order must be smaller')
    p = poly(E,[f(x).diff(x,i)]); c = p.coeffs()[0]
    mat2[i] = -c
    nonzero=False
    for j in range(nvars):
        C = poly(c,[vars[j]]).coeffs()
        if len(C)==2:
            mat1[eq_count,j] = C[0]
            mat2[eq_count,0] += C[0]*vars[j]
            nonzero=True
    if nonzero: eq_count += 1
if eq_count < nvars: raise ValueError('order must be larger')
print((mat1.inv()*mat2))
#check_derivation((mat1.inv()*mat2)[0])
