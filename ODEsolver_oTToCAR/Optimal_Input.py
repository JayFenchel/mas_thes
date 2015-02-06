#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.optimize import fmin
from scipy.optimize import minimize
from time import time
from Modell_kont import OdeModel
from Modell_diskret_new import DiskretOdeModel


class OptimalControlProblem:

    def __init__(self):
        ode_model = OdeModel()
        self.prediction = ode(ode_model.f_t_first).set_integrator('zvode', method='bdf', with_jacobian=False, rtol=1, atol=1)
        self.prediction.set_initial_value(x0, t0).set_jac_params(2.0)

    def cost_function(self, u10):
        u10 = np.array(u10).flat[0]
        us = np.array([u10, 30])  # fmin macht aus u10 ein array, was im weiteren Probleme verursacht
        self.prediction.set_initial_value(x0, t0).set_f_params(us)
        self.prediction.integrate(dt)

        soluti = self.prediction.y.real
        return np.square(soluti.T[2]) + .001*np.square(u10)
        # return soluti.T[2]**2 + .001*u10**2

def cost_function(u1, s):
    us = [u1, 30]
    solution = odeint(myDGL.f, x0, t, args=(us,), rtol=s, atol=s)
    return np.square(solution.T[2][n_sample-1]-0).sum() + .01*np.square(u1)

def cost_function_diskret(u1, xzero):
    us = [u1, 30]
    for i in range(1,5,1):
        solution = np.array(myDiskret.f(xzero, t1/5, us))
        xzero = solution
    return np.square(solution.T[2]-0).sum() + .01*np.square(u1)

def cost_for_t(tend, u, xzero):
    n_sample = 10
    t = np.linspace(0, tend, n_sample)
    solution = odeint(myDGL.f, xzero, t, args=(u,))
    return (solution[-1][1]-0.2)**2


if __name__ == '__main__':
    yaw = 0

    # x0, t0 = np.array([0, 0, 0, yaw, 0, v, 0]), 0
    # t1 = 0.1
    # dt = .1
    # n_sample = (t1-t0)/dt + 1
    # dim_x = x0.shape[0]
    #
    myDGL = OdeModel()
    # myDiskret = DiskretOdeModel()
    # t = np.linspace(t0, t1, n_sample)
    # print(t)
    # myOCP = OptimalControlProblem()
    # u0 = [0.1, 30]
    #
    # time0 = time()
    # uopt = fmin(myOCP.cost_function, u0[0])
    # print(time()-time0)
    # print(uopt)
    #
    # dt = 0.1
    # n_sample = (t1-t0)/dt + 1
    # t = np.linspace(t0, t1, n_sample)
    # time1 = time()
    # uopt = fmin(cost_function, u0[0], args=(.0000000149012,))
    # print(time()-time1)
    # print(uopt)
    #
    # time1 = time()
    # uopt = fmin(cost_function, u0[0], args=(1,))
    # print(time()-time1)
    # print(uopt)
    #
    # time1 = time()
    # uopt = fmin(cost_function_diskret, u0[0], args=(x0,))
    # print(time()-time1)
    # print(uopt)
    #
    # time1 = time()
    # bnds = [(-0.1, 1)]
    # uopt = minimize(cost_function_diskret, u0[0], method='SLSQP',args=(x0,), bounds=bnds)
    # print(uopt)
    # print(time()-time1)
    #
    # t1 = 1
    # dt = 0.01
    # n_sample = (t1-t0)/dt + 1
    # t = np.linspace(t0, t1, n_sample)
    # u = [uopt.x, 30]
    # endsolution = odeint(myDGL.f, x0, t, args=(u,))
    n_sample = 10
    i= 0
    duration = np.zeros_like(np.linspace(1, 4, 10))

    for v in np.linspace(1, 4, 10):
        i
        u = np.array([127*0.0028, 15+20*v])
        x0, t0 = np.array([0, 0, 0, yaw, 0, v, 0]), 0
        tend0 = 0.1
        bnds = [(0, 1)]
        topt = minimize(cost_for_t, tend0, args=(u, x0,), method='SLSQP', bounds=bnds)
        duration[i] = topt.x
        i += 1
        print (topt.x)
        #
        # t = np.linspace(0, topt.x, n_sample)
        # endsolution = odeint(myDGL.f, x0, t, args=(u,))
        #
        # plt.figure(1)
        # plt.subplot(2, 1, 1)
        # plt.title('kontinuierliches System')
        # plt.plot(endsolution.T[0], endsolution.T[1])
        # plt.subplot(2, 1, 2)
        # plt.plot(t, endsolution.T[1])
        # # plt.subplot(3, 1, 3)
        # # plt.plot(t, endsolution.T[5])
        # plt.show()
    v = np.linspace(1, 4, 10)

    plt.plot(v, duration)
    plt.grid()
    plt.show()