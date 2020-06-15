import scipy as sp
import numpy as np
import timeit as tm
from scipy import integrate, special
from scipy.constants import pi
from numpy import exp, inf
from math import sqrt

# Random array of a's, c's and l's
N = 5
values_A = np.random.randint(low=1, high=11, size=[3, N, N])
values_mol = np.random.randint(low=1, high=11, size=[N, 3])
values_r = np.random.randint(low=1, high=11, size=[N, 3])


# Data of sum(den_A) for the system of equations
atom = 0
coeff = []
while atom < N:
    dist = 0
    den_A = []
    while dist < N:
        a_A = values_A[0, atom, dist]
        c_A = values_A[1, atom, dist]
        l_A = values_A[2, atom, dist]
        fact_A = special.factorial2(2*l_A+1, exact=False)
        i_anlt_A = c_A*sqrt((pi/a_A)**3)*fact_A/(2*a_A)**l_A
        den_A.append(i_anlt_A)
        dist += 1
    np.asarray(den_A)
    coeff.append(den_A)
    atom += 1
np.asarray(coeff)


# Data of den_mol for the system of equations
den_mol = []
for file in values_mol:
    a_mol = file[0]
    c_mol = file[1]
    l_mol = file[2]
    fact_mol = special.factorial2(2*l_mol+1, exact=False)
    i_anlt_mol = c_mol*sqrt((pi/a_mol)**3)*fact_mol/(2*a_mol)**l_mol
    den_mol.append(i_anlt_mol)


# Solving the system of linear equations
b = np.linalg.solve(coeff, den_mol)
print("The values of b are:", b)


# Evaluating den_mol_cusp
def evaluating(function):
    function_arr = []
    for file in values_r:
        a_r = file[0]
        c_r = file[1]
        l_r = file[2]
        fact_r = special.factorial2(2*l_r+1, exact=False)
        i_anlt_r = c_r*sqrt((pi/a_r)**3)*fact_r/(2*a_r)**l_r
        function_arr.append(i_anlt_r)
    return print(function, "=", sum(b*function_arr))


evaluating('den_mol_cusp')
