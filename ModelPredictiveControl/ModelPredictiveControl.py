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
v0 = np.eye(T*n, 1)*0

zk = z0
vk = v0
zv_k = np.vstack([zk, vk])

for i in xrange(0, 100):
    for i in xrange(0, T):
        zk[i:i+m], zk[i+m:i+m+n] = u0, x0


    delta_zv, r = QP.solve(zv_k)

    # Schrittweite s in (0,1] bestimmen für die norm(r) minimal ist

    last_r_norm = 10000000000
    for i in np.linspace(1, 0.1, 10):
        zv_help = zv_k + delta_zv
        # TODO Pz < h überprüfen
        r = QP.residual(zv_help)
        r_norm = ((r[:]*r[:]).sum())
        if r_norm < last_r_norm:
            s = i
        else:
            break

    zv_k += s*delta_zv

    print(zv_k[0:m], s, r_norm)

