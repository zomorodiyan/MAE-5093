from sympy import symbols, diff, poly, Function, factorial
import numpy as np
x, h, a, b = symbols('x, h, a, b')
f = Function('f')
def t(j,k): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(x).diff(x, i) for i in range(abs(k)))
eq = t(0.5,4) - 1/8*(3*f(x)+6*t(1,4)-t(2,4))
print(eq)
