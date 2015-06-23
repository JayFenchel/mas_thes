#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
# import cProfile
from ModelPredictiveControl.MyMath import vector_norm
from ModelPredictiveControl.MyMath import backtracking_line_search
from ModelPredictiveControl.MyMath import backtracking_line_search_quick_and_dirty
from ModelPredictiveControl.Systems import SimpleExample
from ModelPredictiveControl.Systems import AirCraft
from ModelPredictiveControl.Systems import qp_from_test
from ModelPredictiveControl.Systems import qp_from_sys
from ModelPredictiveControl.Systems import qp_from_new_sys
from time import time
from ModelPredictiveControl.QuadraticProgram import QuadraticProgram
# import matplotlib.pyplot as plt

# TODO Profiling, mach das mal
# profiler = cProfile.Profile()

# sys = AirCraft()
# QP1 = QuadraticProgram(sys)
QP, AA, bb = qp_from_test()
# exit()
n = QP.n
m = QP.m
T = QP.T  # Planning horizon

# Startwerte AirCraft
x0 = np.array([[0.], [0.], [0.], [400.], [0.]])
u0 = np.array([[0.]])

x0, u0 = QP.x0, QP.u0

z0 = np.zeros([T*(m+n), 1])
for i in range(0, T):
    z0[i*(n+m):i*(m+n)+m], z0[i*(m+n)+m:i*(m+n)+m+n] = u0, x0
v0 = np.ones([T*n, 1])*10
zv_0 = np.vstack([z0, v0])

xk = x0
zv_k = zv_0

print('startwert valide = ', QP.check(xk, zv_k))  # Validität des Startwerts prüfen
zeit = time()
# profiler.enable()
schritte = 25
st_tolerance = 0.0466  # 0.6^2
r_tolerance = 1e-2
x_out = np.zeros([np.shape(xk)[0], schritte])
for schritt in range(schritte):
    x_out[:, schritt:schritt+1] = xk
    QP.kappa = 90
    for zweimal in range(2):
        st, rp_norm, rd_norm = 1, 1, 1
        # Optimize until algorithm fails to go further (st < st_tolerance) or
        # residual is small enough
        while st >= st_tolerance and rp_norm+rd_norm >= r_tolerance:
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
            st = backtracking_line_search(QP.residual_norm, zv_k, delta_zv,
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
            print(zweimal)
        QP.kappa *= 0.1
    # print(zv_k)

    # Ausgabe kostenfunktionswert test-cases
    x_k_plus_eins = zv_k[0:m+n]
    f_von_x = np.dot(QP.g.T, x_k_plus_eins) + .5*np.dot(x_k_plus_eins.T, np.dot(QP.H, x_k_plus_eins))
    print('Value of cost function = %s' % f_von_x)
    x_k_plus_eins[0] += 0.01
    x_k_plus_eins[1] += -0.01
    f_von_x = np.dot(QP.g.T, x_k_plus_eins) + .5*np.dot(x_k_plus_eins.T, np.dot(QP.H, x_k_plus_eins))
    print('Value of cost function geändert = %s' % f_von_x)

    print('x_k_plus1', x_k_plus_eins)
    print('Ax-b', np.dot(AA, x_k_plus_eins) - bb)
    print('Cx - bofx0', np.dot(QP.C, x_k_plus_eins) - QP.b_of_xk(x0))
    print('Adx0- xk+1', np.dot(QP.A, x0) - x_k_plus_eins)
    print(np.dot(QP.g.T, x_k_plus_eins))
    # print(zv_k[0])
    # print(np.dot(sys.B, zv_k[0]))

    # # neues xk berechnen
    # xk = np.dot(QP.A, xk) + QP.B*zv_k[0:m]  #TODO np.dot darf nicht für multiplikation mit skalaren genommen werden
    # # z_k shiften
    # zv_k[0:(n+m)*(T-1)] = zv_k[n+m:(n+m)*T]
    # # neues z_k[T] hinten anhängen, u[T] ist nicht ganz korrekt, aber kein
    # # u[T+1] vorhanden
    # zv_k[(n+m)*(T-1)+m:(n+m)*T] = np.dot(QP.A, zv_k[(n+m)*(T-1)+m:(n+m)*T])\
    #                     + QP.B*zv_k[(T-1)*(m+n):(T-1)*(m+n)+m]
    # v_k shiften
    zv_k[(n+m)*T:(n+m)*T+n*(T-1)] = zv_k[(n+m)*T+n:(n+m)*T+n*T]
    # neues v_k[T] hinten anhängen
    zv_k[(n+m)*T+n*(T-1):(n+m)*T+n*(T)] = np.ones([n, 1])*100
    print('next xk', xk)
print(time()-zeit)
# profiler.disable()
# profiler.dump_stats('profile')
print(x_out)