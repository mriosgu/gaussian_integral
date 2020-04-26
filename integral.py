import scipy as sp
import numpy as np
import timeit as tm
from scipy import integrate, special
from scipy.constants import pi
from numpy import exp, inf
from math import sqrt

# Random array of a's, c's and l's
values = np.random.randint(low=1, high=11, size=3)
a = values[0]
c = values[1]
l = values[2]

# Numerical integration
def f(r, a, c, l):
    return c*r**(2*l)*exp(-a*r**2)*4*pi*r**2
i_num, err = integrate.quad(f, 0, inf, args=(a, c, l))
print("Integral_num = {:f} (+-{:g})".format(i_num, err))

# Analytical integration
fact = special.factorial2(2*l+1, exact=False)
i_anlt = c*sqrt((pi/a)**3)*fact/(2*a)**l
print("Integral_anlt = ", i_anlt)

# Test if numerical and analytical integration yields the same result
if np.allclose(i_num, i_anlt) == True:
    print("Same result")
else:
    print("Something is wrong")

# Test execution time 
print("time_num = ", tm.timeit('i_num', number=10000, globals=globals()))
print("time_anlt = ", tm.timeit('i_anlt', number=10000, globals=globals()))
