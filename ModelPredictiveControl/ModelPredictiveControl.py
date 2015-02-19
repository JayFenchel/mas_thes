#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from MyMath import vector_norm
from Systems import SimpleExample
from Systems import AirCraft
from QuadraticProgram import QuadraticProgram
import matplotlib.pyplot as plt


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
for i in xrange(0, T):
    zk[i*(n+m):i*(m+n)+m], zk[i*(m+n)+m:i*(m+n)+m+n] = u0, x0
zk[0] = 0
zv_k0 = np.vstack([zk, vk])
xk = x0
zv_k = zv_k0
print 'startwert valide = ',QP.check(zv_k)  # Validität des Startwerts prüfen
for i in xrange(0, 100):
    delta_zv = QP.solve(xk, zv_k)

    # Schrittweite s in (0,1] bestimmen für die norm(r) minimal ist
    # backtracking line search
    f_x = np.square(np.vstack(QP.residual(zv_k))).sum()
    f_xp = np.zeros([100, 1])
    for i in range (0, 100, 1):
        f_xp[i] = np.square(np.vstack(QP.residual(zv_k + (100.-i)*delta_zv/100.))).sum()
        # if not QP.check(zv_k + (100.-i)*delta_zv/100.):
        #     f_xp[i] = 0
        # print ((100.-i)/100., QP.check(zv_k + (100-i)*delta_zv/100))
    plt.plot(np.linspace(1, 0, 100), f_xp)
    plt.grid()
    # plt.show()
    testschritt = delta_zv * .0000001
    DELTA = np.eye(np.shape(zv_k)[0])*0.000001
    HILF = (zv_k + DELTA)
    # print(HILF.T[0:1].T)
    delta_f = np.zeros([np.shape(zv_k)[0], 1])
    for k in range(0, np.shape(zv_k)[0]):
        # print (np.square(np.vstack(QP.residual(HILF.T[0:1].T))).sum() - f_x)/np.square(0.00001)
        delta_f[k] = (np.square(np.vstack(QP.residual(HILF.T[k:k+1].T))).sum() - f_x)/0.000001
    nabla_f = (np.square(np.vstack(QP.residual(zv_k + testschritt))).sum() - f_x)/np.sqrt(np.square(testschritt).sum())
    alpha = 0.4
    beta = 0.6
    st = 1
    print np.dot(delta_f.T, delta_zv)
    while np.square(np.vstack(QP.residual(zv_k + st*delta_zv))).sum() > f_x + alpha*st*np.dot(delta_f.T, delta_zv):
        st = beta*st
        print st
    if QP.check(zv_k + st*delta_zv):
        print('Valid step possible')
        zv_k += st*delta_zv
    else:
        st = 0
    # print(zv_k)
    for i in range (0,T):
        print(zv_k[i*(m+n):i*(m+n)+m])
    rd, rp = QP.residual(zv_k)
    rd_norm = np.square(rd[:]).sum()
    rp_norm = np.square(rp[:]).sum()
    r_norm = rd_norm + rp_norm
    print(st, rp_norm, rd_norm)
print zv_k
# print rp
# xk = x0
# zv_k = zv_k0
#
# for i in xrange(0, 20):
#
#
#     delta_zv, r = QP.solve_own(xk, zv_k)
#
#     # Schrittweite s in (0,1] bestimmen für die norm(r) minimal ist
#
#     last_r_norm = 10000000000
#     s = 0
#     for i in np.linspace(1, .1, 10):
#         zv_help = zv_k + i*delta_zv
#         if QP.check(zv_help):
#             r = QP.residual(zv_help)
#             r_norm = ((r[:]*r[:]).sum())
#             if r_norm < last_r_norm:
#                 s = i
#                 last_r_norm = r_norm
#             else:
#                 break
#     if s == 0:
#         print('No valid step possible')
#
#     zv_k += s*delta_zv
# for i in range (0,T):
#     print(zv_k[i:(i+1)*m])
# print(s, r_norm)
#
