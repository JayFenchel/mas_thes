#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from MyMath import matrix_diag
from MyMath import solve_lin_gs



class QuadraticProgram:

    def __init__(self, sys):

        n = sys.n
        m = sys.m
        T = sys.T  # Planning horizon, Anzahl der Schritte

        self.delta_t = sys.delta_t  # Länge der Schritte # TODO richtige Zeitschitte einbauen

        self.n = n
        self.m = m
        self.T = T
        self.A = sys.A
        self.S = sys.S
        self.Fx = sys.Fx

        # Cost function
        H = np.eye(T*(sys.n+sys.m), T*(sys.n+sys.m))
        H[0:m].T[0:m] = sys.R
        QSR = np.hstack([np.vstack([sys.Q, sys.S.T]), np.vstack([sys.S, sys.R])])
        for i in range(1, T):
            # TODO wegen dem Transponieren auf Symmetrie achten
            H[m+(i-1)*(m+n):m+i*(m+n)].T[m+(i-1)*(m+n):m+i*(m+n)] = QSR
        H[m+(T-1)*(m+n):m+(T-1)*(m+n)+n].T[m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = sys.Qf
        self.H = H

        g = np.eye(T*(m+n), 1)
        for i in range(0,T):
            g[i*(n+m):(i+1)*(n+m)] = np.vstack([sys.r, sys.q])
        g[(T-1)*(n+m)+m:T*(n+m)] = sys.qf
        self.g = g

        # Equality constraints
        C = np.eye(T*n, T*(n+m))*0
        C[0:n].T[0:m+n] = np.vstack([-sys.B.T, np.eye(n, n)])
        for i in range(1, T):
            C[i*n:(i+1)*n].T[m+(i-1)*(m+n):m+i*(m+n)+n] = np.vstack([-sys.A.T, -sys.B.T, np.eye(n, n)])
        self.C = C

        self.b = np.eye(T*n, 1)

        # Inequality constraints
        P = np.eye(T*2*(n+m)+n, T*(n+m))*0
        P[0:2*(n+m)].T[0:m] = sys.Fu.T
        for i in range(1, T):
            P[i*2*(n+m):(i+1)*2*(n+m)].T[m+(i-1)*(m+n):m+i*(m+n)] = np.vstack([sys.Fx.T, sys.Fu.T])
        P[T*2*(n+m):T*2*(n+m)+n].T[m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = sys.Ff.T
        self.P = P

        h = np.eye(2*T*(m+n)+n, 1)
        for i in range(0,T):
            h[i*2*(m+n):(i+1)*2*(m+n)] = sys.f
        h[2*T*(m+n):2*T*(m+n)+n] = sys.ff
        self.h = h
        self.v0 = np.eye(T*n, 1)*0



    def solve(self, zv_k):

        T = self.T
        n = self.n
        m = self.m

        xk = zv_k[m:m+n]
        self.b[0:n] = np.dot(self.A, xk)
        self.g[0:m] += 2*np.dot(self.S.T, xk)
        self.h[0:2*(n+m)] += -np.dot(self.Fx, xk)

        self.kappa = 10  # >0 barrier parameter

        self.d = np.eye(2*self.T*(self.m+self.n)+n, 1)
        self.d[:] = 1/(self.h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))



        Phi = 2*self.H + self.kappa*np.dot(np.dot(self.P.T, matrix_diag(self.d)), self.P)

        r = self.residual(zv_k)
        SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.eye(self.C.shape[0], self.C.shape[0])*0])])


        lsg = np.linalg.solve(SS, -r)
        return lsg, r

    def solve_own(self, zv_k):

        T = self.T
        n = self.n
        m = self.m

        xk = zv_k[m:m+n]
        self.b[0:n] = np.dot(self.A, xk)
        self.g[0:m] += 2*np.dot(self.S.T, xk)
        self.h[0:2*(n+m)] += -np.dot(self.Fx, xk)

        self.kappa = 10  # >0 barrier parameter

        self.d = np.eye(2*self.T*(self.m+self.n)+n, 1)
        self.d[:] = 1/(self.h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))



        Phi = 2*self.H + self.kappa*np.dot(np.dot(self.P.T, matrix_diag(self.d)), self.P)

        r = self.residual(zv_k)
        SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.eye(self.C.shape[0], self.C.shape[0])*0])])


        lsg = np.linalg.solve(SS, -r)
        lsg = solve_lin_gs(SS, -r)
        return lsg, r

    def residual(self, zv_k):
        rd = 2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, self.d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):])
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - self.b

        return np.vstack([rd, rp])

    def check(self, zv_k):
        if((np.dot(self.P, zv_k[0:self.T*(self.m+self.n)]) - self.h) > 0).all():
            return False
        return True



