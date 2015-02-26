#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from MyMath import matrix_diag
from MyMath import solve_lin_gs
from MyMath import solve_lin_gs_structured



class QuadraticProgram:


    def set_sys_dynamics(self, A, B):

        # Equality constraints
        C = np.zeros([T*n, T*(n+m)])
        C[0:n, 0:m+n] = np.hstack([-sys.B, np.eye(n, n)])
        for i in range(1, T):
            C[i*n:(i+1)*n, m+(i-1)*(m+n):m+i*(m+n)+n] = np.hstack([-sys.A, -sys.B, np.eye(n, n)])
        self.C = C

        self.b = np.zeros([T*n, 1])
        pass

    def set_weighting(self):
        pass

    def set_constraints(self, Fu, fu, Fx, fx, Ff, ff):

        if self.n is None:
            print('You have to set_sys_dynamics() fist.')
            exit()
        T, n, m = 10, 5, 1

        # Inequality constraints
        n_Fu = np.shape(Fu)[0]
        P = np.zeros([T*n_Fu+np.shape(Ff)[0], T*(n+m)])
        P[0:n_Fu, 0:m] = Fu
        for i in range(1, T):
            Hilf = np.hstack([Fx, Fu])
            P[i*n_Fu:(i+1)*n_Fu, m+(i-1)*(m+n):m+i*(m+n)] = Hilf

        P[T*n_Fu:T*n_Fu+np.shape(Ff)[0], m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = Ff
        self.P = P

        f = np.vstack([fx, fu]) # stacken - anderer branch
        f = fx + fu  # TODO fraglich

        h = np.zeros([T*np.shape(f)[0]+np.shape(ff)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f
        h[T*np.shape(f)[0]:T*np.shape(f)[0]+np.shape(ff)[0]] = ff
        self.h = h

    def __init__(self):
        self.n = None

    # def __init__(self, sys):
    #
    #
    #     self.kappa = 90#>0 barrier parameter
    #
    #     n = sys.n
    #     m = sys.m
    #     T = sys.T  # Planning horizon, Anzahl der Schritte
    #
    #     self.delta_t = sys.delta_t  # Länge der Schritte # TODO richtige Zeitschitte einbauen
    #
    #     self.n = n
    #     self.m = m
    #     self.T = T
    #     self.A = sys.A
    #     self.B = sys.B
    #     self.r = sys.r
    #     self.f = sys.f
    #     self.S = sys.S
    #     self.Fx = sys.Fx
    #
    #     # Cost function
    #     H = np.eye(T*(n+m), T*(n+m))
    #     H[0:m, 0:m] = sys.R
    #     QSR = np.hstack([np.vstack([sys.Q, sys.S.T]), np.vstack([sys.S, sys.R])])
    #     for i in range(1, T):
    #         H[m+(i-1)*(m+n):m+i*(m+n), m+(i-1)*(m+n):m+i*(m+n)] = QSR
    #     H[m+(T-1)*(m+n):m+(T-1)*(m+n)+n, m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = sys.Qf
    #     self.H = H
    #
    #     g = np.zeros([T*(m+n), 1])
    #     for i in range(0, T):
    #         g[i*(n+m):(i+1)*(n+m)] = np.vstack([sys.r, sys.q])
    #     g[(T-1)*(n+m)+m:T*(n+m)] = sys.qf
    #     self.g = g
    #
    #     # Equality constraints
    #     C = np.zeros([T*n, T*(n+m)])
    #     C[0:n, 0:m+n] = np.hstack([-sys.B, np.eye(n, n)])
    #     for i in range(1, T):
    #         C[i*n:(i+1)*n, m+(i-1)*(m+n):m+i*(m+n)+n] = np.hstack([-sys.A, -sys.B, np.eye(n, n)])
    #     self.C = C
    #
    #     self.b = np.zeros([T*n, 1])

    def solve(self, xk, zv_k):

        T = self.T
        n = self.n
        m = self.m

        self.b[0:self.n] = np.dot(self.A, xk)
        self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)
        self.h[0:np.shape(self.Fx)[0]] = self.f - np.dot(self.Fx, xk)

        d = np.zeros([np.shape(self.P)[0], 1])
        d[:] = 1/(self.h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))

        Phi = 2*self.H + self.kappa*np.dot(np.dot(self.P.T, matrix_diag(d**2)), self.P)  # TODO vernünftiges Quadrieren

        # print(m)
        # print(self.P[0:m+n+7].T[0:m+n+3]).T
        # print(np.linalg.inv(Phi))

        r = np.vstack(self.residual(xk, zv_k))
        SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.eye(self.C.shape[0], self.C.shape[0])*0])])

        # v = zv_k[self.T*(self.m+self.n):]
        # lsg = solve_lin_gs_structured(Phi, r, self.A, self.B, self.C, T, m, n, v)
        # print lsg[100:]

        lsg = np.linalg.solve(SS, -r)

        # print abs(-r - np.dot(SS, lsg)).sum()
        # print lsg[100:]
        return lsg
    #
    # def solve_own(self, xk, zv_k):
    #
    #     T = self.T
    #     n = self.n
    #     m = self.m
    #
    #     self.b[0:n] = np.dot(self.A, xk)
    #     self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)
    #     self.h[0:np.shape(self.Fx)[0]] = self.f - np.dot(self.Fx, xk)
    #
    #     self.kappa = .001  # >0 barrier parameter
    #
    #     self.d = np.eye(np.shape(self.P)[0], 1)
    #     self.d[:] = 1/(self.h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))
    #
    #
    #
    #     Phi = 2*self.H + self.kappa*np.dot(np.dot(self.P.T, matrix_diag(self.d)), self.P)
    #
    #     r = self.residual(zv_k)
    #     SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.eye(self.C.shape[0], self.C.shape[0])*0])])
    #
    #     lsg = solve_lin_gs(SS, -r)
    #     return lsg, r

    def residual(self, xk, zv_k):
        h=self.h
        h[0:np.shape(self.Fx)[0]] = self.f - np.dot(self.Fx, xk)

        d = np.zeros([np.shape(self.P)[0], 1])
        d[:] = 1/(h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))

        # print 'term12',np.square(2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, d)).sum()
        rd = 2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):])
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - self.b

        if not self.check(zv_k):
            return rd + 100000000000000000000000000000000000000, rp + 100000000000000000000000000000000000000

        return rd, rp

    def check(self, zv_k):
        return ((np.dot(self.P, zv_k[0:self.T*(self.m+self.n)]) - self.h) < 0).all()



