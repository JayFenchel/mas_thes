#!/usr/bin/env python2
import ctypes
import time
import numpy as np
from OdeModel import OdeModel
import timeit

from numpy.ctypeslib import ndpointer
mylib = ctypes.cdll.LoadLibrary("hello.so")
mylib.func.restype = ndpointer(dtype=ctypes.c_float, shape=(7,))

model = OdeModel()


dummy1=time.clock()
for i in xrange(10000):
    model.f_discrete(np.array([0,0,0,0.5,0,2,0]),0.06,np.array([0.4,45]))
print time.clock()-dummy1


dummy1=time.clock()
for i in xrange(10000):
    x = (ctypes.c_float*7)()
    x[0]=0
    x[1]=0
    x[2]=0
    x[3]=0.5
    x[4]=0
    x[5]=2
    x[6]=0
    
    u = (ctypes.c_float*2)()
    
    u[0]=0.4
    u[1]=45
    
    xp=ctypes.cast(x, ctypes.POINTER(ctypes.c_float))
    up=ctypes.cast(u, ctypes.POINTER(ctypes.c_float))
    mylib.func(xp,up,ctypes.c_float(0.06))
print time.clock()-dummy1





