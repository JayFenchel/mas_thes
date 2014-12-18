#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:00:24 2014

@author: jayf
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.integrate import odeint

from time import sleep, time

x0, t0 = np.array([0, 0, 0, 0.1, 1]), 0
t1 = 1
dt = .01
n_sample = (t1-t0)/dt + 1
dim_x = x0.shape[0]

def f(x, t, p):
    
    print(x)
    
    dx1 = x[4]*np.cos(x[2])
    dx2 = x[4]*np.sin(x[2])
    dx3 = x[3]
    dx4 = 0
    dx6 = 0
## TODO Meine Modellgleichungen aus Matlab verwenden

#    dx0 = x[5]*np.cos(x[2]+x[4]);
#    dx1 = x[5]*np.sin(x[2]+x[4]);
#    dx2 = x[3];
#    dx3 = 1 / Jz * (Fvq * lv * cos(x[6]) + Fvl * lv * sin(x[6]) - Fhq * lh);
#    dx4 = (1 / (m * x[5])) * (Fvq * cos(x[6]-x[4]) + Fvl * sin (x[6]-x[4]) + Fhq * cos(x[4]) + FK * sin(x[4]) - Fhl * sin(x[4])) - x[3];
#    dx5 = (1/m) * (Fvl * cos(x[6] - x[4]) + Fhq * sin(x[4]) + Fhl * cos(x[4]) - Fvq * sin(x[6]-x[4]) - FR - FL - FK * cos(x[4]));
#    dx6 = -x[6]/T + u0/T;
#    return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]
    return [dx1, dx2, dx3, dx4, dx6]
    
time0 = time()   
r = ode(f).set_integrator('zvode', method='bdf', with_jacobian=False)


#~ for i in range(50):  
    #~ r.set_initial_value(x0, t0).set_f_params(1.0).set_jac_params(2.0)
  ###  xlsg=np.matrix(x0)
  ###  tlsg=np.matrix(t0)
    #~ xlsg=np.zeros((n_sample, dim_x))
    #~ tlsg=np.zeros(n_sample)
    #~ xlsg[0] = x0
    #~ tlsg[0] = t0
    #~ 
    #~ k = 1
    #~ while r.successful() and r.t < t1:
        #~ r.integrate(r.t+dt)
      ###  xlsg = np.vstack([xlsg, r.y.real])
      ###  tlsg = np.vstack([tlsg, r.t])
        #~ xlsg[k] = r.y.real
        #~ tlsg[k] = r.t
        #~ 
        #~ k += 1
    #~ xlsg = xlsg.T
#~ 
#~ print(time()-time0)

for i in range(1):  
    t = np.linspace(0, .1, 100)
    y0=x0
    soln = odeint(f, y0, t, args = (0,))
#~ 
#~ print(xlsg)
#~ print(tlsg)
# plt.figure(1)
# plt.subplot(3,1,1)
# plt.title('kontinuierliches System')
# plt.plot(tlsg,xlsg[0])
# plt.subplot(3,1,2)
# plt.plot(tlsg,xlsg[1])
# plt.subplot(3,1,3)
# plt.plot(tlsg,xlsg[2])
# plt.show()
input()
