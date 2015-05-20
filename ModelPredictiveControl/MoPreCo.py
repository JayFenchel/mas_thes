#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
# import cProfile
from ModelPredictiveControl.MyMath import vector_norm
from ModelPredictiveControl.MyMath import backtracking_line_search
from ModelPredictiveControl.MyMath import backtracking_line_search_better
from ModelPredictiveControl.MyMath import backtracking_line_search_quick_and_dirty
from ModelPredictiveControl.Systems import SimpleExample
from ModelPredictiveControl.Systems import AirCraft
from ModelPredictiveControl.Systems import qp_from_sys
from ModelPredictiveControl.Systems import qp_from_new_sys
from time import time
from ModelPredictiveControl.QuadraticProgram import QuadraticProgram
# import matplotlib.pyplot as plt

# TODO Profiling, mach das mal
# profiler = cProfile.Profile()

# sys = AirCraft()
# QP1 = QuadraticProgram(sys)
QP = qp_from_new_sys()

n = QP.n
m = QP.m
T = QP.T  # Planning horizon

# Startwerte AirCraft
x0 = np.array([[0.], [0.], [0.], [400.], [0.]])
x0 = np.array([[0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [400.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.]])
u0 = np.array([[0.]])
z0 = np.eye(T*(n+m), 1)*0
v0 = np.eye(T*n, 1)*0

# # Startwere SimpleSys
# x0 = np.array([[0], [1], [0]])
# u0 = np.array([[0],[0]])
# z0 = np.eye(T*(n+m), 1)*0
# v0 = np.eye(T*n, 1)*0

zk = z0
vk = v0
for i in range(0, T):
    zk[i*(n+m):i*(m+n)+m], zk[i*(m+n)+m:i*(m+n)+m+n] = u0, x0
zk[0] = 0
zv_k0 = np.vstack([zk, vk])
xk = x0
zv_k = zv_k0
# zv_k = np.array([[  2.09719785e-01],
#                  [ -2.58929379e-01], [ -3.01636240e-01], [ -9.40136999e-01], [  3.99622527e+02], [  2.09719785e-01],
#                  [ -7.75719981e-02],
#                  [ -1.34373473e-01], [ -3.47773039e-01], [  4.73557241e-01], [  3.90842471e+02], [ -7.75719981e-02],
#                  [  4.54511936e-02],
#                  [ -3.72957122e-03], [ -2.35142989e-01], [ -4.95019233e-03], [  3.76180207e+02], [  4.54511936e-02],
#                  [ -1.12575574e-03],
#                  [ -3.89703411e-04], [ -2.33473482e-01], [  8.05178514e-03], [  3.61280199e+02], [ -1.12575574e-03],
#                  [  5.67869071e-04],
#                  [  6.44327917e-04], [ -2.31969075e-01], [ -1.04174820e-03], [  3.46351916e+02], [  5.67869071e-04],
#                  [ -6.50539423e-04],
#                  [  7.71624133e-04], [ -2.31554793e-01], [  2.13365617e-03], [  3.31449839e+02], [ -6.50539423e-04],
#                  [ -3.69509271e-04],
#                  [  1.02268345e-03], [ -2.30733896e-01], [  1.18874870e-03], [  3.16575201e+02], [ -3.69509271e-04],
#                  [ -2.35626886e-03],
#                  [  3.36700235e-03], [ -2.27404301e-01], [  9.71526302e-03], [  3.01742206e+02], [ -2.35626886e-03],
#                  [ -8.74586923e-03],
#                  [  1.33422283e-02], [ -2.13452843e-01], [  3.72216289e-02], [  2.87039876e+02], [ -8.74586923e-03],
#                  [ -6.45648477e-03],
#                  [  1.78250506e-02], [ -1.99073137e-01], [  2.09018996e-02], [  2.72805566e+02], [ -6.45648477e-03],
#                  # v0
#                  [  1.11578458e+05], [ -1.25767521e+05], [  1.01824826e+04], [ -1.81812211e+03], [ -4.21509254e+02],
#                  [  6.37119096e+04], [ -2.53520706e+04], [ -8.67608160e+03], [ -1.59105832e+03], [  1.44231753e+02],
#                  [  4.90375252e+02], [ -6.77530761e+02], [  5.91885827e+02], [ -1.37048300e+03], [ -6.85235834e+01],
#                  [ -4.92505765e+04], [  4.50713939e+04], [ -4.23958379e+02], [ -1.16074309e+03], [  1.98777238e+00],
#                  [ -2.93838565e+04], [  2.68017679e+04], [ -1.20689973e+02], [ -9.62014279e+02], [ -1.20135709e+00],
#                  [ -2.52785117e+04], [  2.31837881e+04], [ -1.64878625e+02], [ -7.74317474e+02], [  7.67316948e-01],
#                  [ -1.87729473e+04], [  1.72649510e+04], [ -1.29120232e+02], [ -5.97633303e+02], [ -7.36750118e-01],
#                  [ -1.26381971e+04], [  1.18373201e+04], [ -1.43266652e+02], [ -4.31941490e+02], [ -1.26801443e+00],
#                  [ -6.04973336e+03], [  6.49401320e+03], [ -3.01319024e+02], [ -2.77211260e+02], [  9.50242423e+00],
#                  [ -4.21055884e+03], [  4.53808858e+03], [ -2.37228200e+02], [ -1.33346052e+02], [  2.02084617e+01]])
print('startwert valide = ', QP.check(xk, zv_k))  # Validität des Startwerts prüfen
zeit = time()
# profiler.enable()
for schritt in range(15):
    for i in range(0, 6):
        zeits = time()
        delta_zv = QP.solve(xk, zv_k)
        print('solve', 5*(time()-zeits))
        zeit1=time()
        # st = backtracking_line_search(QP.residual_norm, zv_k, delta_zv,
        #                               args=(xk, ))
        # print('first', time()-zeit1)
        # zeit2 = time()
        # st = backtracking_line_search_quick_and_dirty(QP.residual_norm, zv_k, delta_zv,
        #                               args=(xk, ))
        # print('quick and dirty', time()-zeit2)
        # zeit2 = time()
        st = backtracking_line_search_better(QP.residual_norm, zv_k, delta_zv,
                                      args=(xk, ))
        # print('better', time()-zeit2)
        if QP.check(xk, zv_k + st*delta_zv):
            print('Valid step possible')
            zv_k += st*delta_zv
        else:
            st = 0
        # print(zv_k)
        for i in range (0,T):
            print(zv_k[i*(m+n):i*(m+n)+m])
        rd, rp = QP.residual(xk, zv_k)
        rd_norm = np.square(rd[:]).sum()
        rp_norm = np.square(rp[:]).sum()
        r_norm = rd_norm + rp_norm
        print(st, rp_norm, rd_norm)
    # print(zv_k)
    # print(zv_k[0])
    # print(np.dot(sys.B, zv_k[0]))
    xk, zv_k[0:(n+m)*(T-1)] = np.dot(sys.A, xk) + sys.B*zv_k[0], zv_k[n+m:(n+m)*T]  #TODO np.dot darf nicht für multiplikation mit skalaren genommen werden
    print('next xk', xk)
print(time()-zeit)
# profiler.disable()
# profiler.dump_stats('profile')