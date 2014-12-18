#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 16:00:24 2014

@author: jayf
"""
import numpy as np
#~ import matplotlib.pyplot as plt
from scipy.integrate import ode

from time import time



def f(x):
    delta_t = 0.1
    x1kplus1 = x[0] + (x[4]*np.cos(x[2])) * delta_t
    x2kplus1 = x[1] + (x[4]*np.sin(x[2])) * delta_t
    x3kplus1 = x[2] + (x[3]) * delta_t
    x4kplus1 = x[3] + (0) * delta_t
    x6kplus1 = x[4] + (0) * delta_t
    
#    dx0 = x[5]*np.cos(x[2]+x[4]);
#    dx1 = x[5]*np.sin(x[2]+x[4]);
#    dx2 = x[3];
#    dx3 = 1 / Jz * (Fvq * lv * cos(x[6]) + Fvl * lv * sin(x[6]) - Fhq * lh);
#    dx4 = (1 / (m * x[5])) * (Fvq * cos(x[6]-x[4]) + Fvl * sin (x[6]-x[4]) + Fhq * cos(x[4]) + FK * sin(x[4]) - Fhl * sin(x[4])) - x[3];
#    dx5 = (1/m) * (Fvl * cos(x[6] - x[4]) + Fhq * sin(x[4]) + Fhl * cos(x[4]) - Fvq * sin(x[6]-x[4]) - FR - FL - FK * cos(x[4]));
#    dx6 = -x[6]/T + u0/T;
#    return [dx0, dx1, dx2, dx3, dx4, dx5, dx6]
    return [x1kplus1, x2kplus1, x3kplus1, x4kplus1, x6kplus1]
    
time0 = time()   

t1 = 1
dt = .01
x0, t0 = [0, 0, 0, 0.1, 1], 0

for i in range(50):  
    xlsg=np.zeros((101, 5))
    tlsg=np.zeros(101)
    xlsg[0] = x0
    tlsg[0] = t0
    
    k = 1
    while (t0 + k*dt) <= t1:
        xlsg[k] = f(xlsg[k-1])
        tlsg[k] = t0 + k*dt
        k += 1
#        xlsg = np.vstack([xlsg, x])
#        tlsg = np.vstack([tlsg, t])
#    xlsg = xlsg.T

print(time()-time0)
#plt.figure(1)
#plt.subplot(3,1,1)
#plt.title('diskretes System')
#plt.plot(tlsg,xlsg[0])
#plt.subplot(3,1,2)
#plt.plot(tlsg,xlsg[1])
#plt.subplot(3,1,3)
#plt.plot(tlsg,xlsg[2])
#plt.show()
input()
