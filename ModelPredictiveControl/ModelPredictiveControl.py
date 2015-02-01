#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from MyMath import vector_norm
from Systems import SimpleExample
from QuadraticProgram import QuadraticProgram



sys = SimpleExample()
QP = QuadraticProgram(sys)

n = sys.n
m = sys.m
T = sys.T  # Planning horizon


x0 = np.array([[0], [1], [0]])
u0 = np.array([[0],[0]])
z0 = np.eye(T*(n+m), 1)*0
for i in xrange(0, T):
    z0[i:i+m], z0[i+m:i+m+n] = u0, x0
Phi = QP.solve(z0)




print(Phi)

