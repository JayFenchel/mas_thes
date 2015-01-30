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
QP.solve(x0)



# print(QP.C)

