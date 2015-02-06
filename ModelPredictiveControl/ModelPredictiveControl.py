#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from MyMath import vector_norm
from Systems import SimpleExample
from Systems import AirCraft
from QuadraticProgram import QuadraticProgram



sys = AirCraft()
QP = QuadraticProgram(sys)

n = sys.n
m = sys.m
T = sys.T  # Planning horizon

# Startwerte AirCraft
x0 = np.array([[0], [0], [0], [400], [0]])
u0 = np.array([[0]])
z0 = np.eye(T*(n+m), 1)*0
v0 = np.eye(T*n, 1)*0

# # Startwere SimpleSys
# x0 = np.array([[0], [1], [0]])
# u0 = np.array([[0],[0]])
# z0 = np.eye(T*(n+m), 1)*0
# v0 = np.eye(T*n, 1)*0

zk = z0
vk = v0
zv_k = np.vstack([zk, vk])

for i in xrange(0, 20):
    for i in xrange(0, T):
        zk[i:i+m], zk[i+m:i+m+n] = u0, x0


    delta_zv, r = QP.solve(zv_k)

    # Schrittweite s in (0,1] bestimmen für die norm(r) minimal ist

    last_r_norm = 10000000000
    for i in np.linspace(1, .1, 10):
        zv_help = zv_k + i*delta_zv
        print(QP.check(zv_help))
        r = QP.residual(zv_help)
        r_norm = ((r[:]*r[:]).sum())
        if r_norm < last_r_norm:
            s = i
            last_r_norm = r_norm
        else:
            break

    zv_k += s*delta_zv
for i in range (0,T):
    print(zv_k[i:(i+1)*m])
print(s, r_norm)

for i in xrange(0, 20):
    for i in xrange(0, T):
        zk[i:i+m], zk[i+m:i+m+n] = u0, x0


    delta_zv, r = QP.solve_own(zv_k)

    # Schrittweite s in (0,1] bestimmen für die norm(r) minimal ist

    last_r_norm = 10000000000
    for i in np.linspace(1, .1, 10):
        zv_help = zv_k + i*delta_zv
        # TODO Pz < h überprüfen
        print(QP.check(zv_help))
        r = QP.residual(zv_help)
        r_norm = ((r[:]*r[:]).sum())
        if r_norm < last_r_norm:
            s = i
            last_r_norm = r_norm
        else:
            break

    zv_k += s*delta_zv
for i in range(0, T):
    print(zv_k[i:(i+1)*m])
print(s, r_norm)

