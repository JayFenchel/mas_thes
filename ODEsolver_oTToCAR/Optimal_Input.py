#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.optimize import fmin
from time import time
from Modell_kont import OdeModel


class OptimalControlProblem:

    def __init__(self):
        ode_model = OdeModel()
        self.prediction = ode(ode_model.f_t_first).set_integrator('zvode', method='bdf', with_jacobian=False)

    def cost_function(self, u10):
        u10 = np.array(u10).flat[0]
        us = np.array([u10, 30])  # fmin macht aus u10 ein array, was im weiteren Probleme verursacht
        self.prediction.set_initial_value(x0, t0).set_f_params(us).set_jac_params(2.0)
        self.prediction.integrate(dt)
        soluti = self.prediction.y.real
        return np.square(soluti.T[2])


def cost_function(u1):
    us = [u1, 30]
    solution = odeint(myDGL.f, x0, t, args=(us,))

    return np.square(solution.T[2][1:n_sample]).sum()

x0, t0 = np.array([0, 0, 0.2, 0, 0, 1, 0]), 0
t1 = 1
dt = .5
n_sample = (t1-t0)/dt + 1
dim_x = x0.shape[0]

myDGL = OdeModel()
t = np.linspace(0, t1, 21)
myOCP = OptimalControlProblem()
time0 = time()

u0 = np.array([0.1, 30])
uopt = fmin(myOCP.cost_function, u0[0])
print(uopt)

print(time()-time0)

u = [uopt, 30]
solution = odeint(myDGL.f, x0, t, args=(u,))


plt.figure(1)
plt.subplot(3, 1, 1)
plt.title('kontinuierliches System')
plt.plot(solution.T[0], solution.T[1])
plt.subplot(3, 1, 2)
plt.plot(t, solution.T[2])
plt.subplot(3, 1, 3)
plt.plot(t, solution.T[5])
plt.show()