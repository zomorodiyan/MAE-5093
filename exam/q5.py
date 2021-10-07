from sympy import symbols, Function, factorial, diff, poly
import matplotlib.pyplot as plt
import numpy as np
y, t, h = symbols('y, t, h')
f = Function('f')
def t(j,k,symb): #taylor f(x+jh) to the kth expression
    return sum((h*j)**i/factorial(i)*f(symb).diff(symb, i) for i in range(abs(k)))
k1 = h*f(x)
