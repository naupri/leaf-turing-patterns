# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:56:16 2026

@author: dgten
"""

import pde
#from sympy import symbols

grid = pde.CartesianGrid([[0, 30], [0, 6]], [30,6], periodic = [True, True])

#k1, k2, k3, k4, k5, k6, k7, k8, k9 = symbols('k1 k2 k3 k4 k5 k6 k7 k8, k9')
k1 = 0
k2 = -0.5
k3 = -1
k4 = -0.5
k5 = -1
k6 = -1
k7 = 1
k8 = 0
k9 = 0 
 #[0, -0.5, -1, -0.5, -1, -1, 1, 0, 0]
Dm = 10
DA = 40

HDZIPIII = pde.ScalarField.random_normal(grid, mean=0)
miRKANARF = pde.ScalarField.random_normal(grid, mean=0.1)
AStasiARF = pde.ScalarField.random_normal(grid, mean=0)

state = pde.FieldCollection([HDZIPIII, miRKANARF, AStasiARF])

eq = pde.PDE({"HDZIPIII": "0 * HDZIPIII + -0.5 * miRKANARF + -1 * AStasiARF",
              "miRKANARF" : "-0.5 * HDZIPIII + -1 * miRKANARF + -1 * AStasiARF + 7.5 * laplace(miRKANARF)",
             "AStasiARF" : "1 * HDZIPIII + 0 * miRKANARF + 0 * AStasiARF + 75 * laplace(AStasiARF)"})

result = eq.solve(state, t_range=100, dt=1)
result.plot()