#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from ModelPredictiveControl.MyMath import matrix_diag
from ModelPredictiveControl.MyMath import solve_lin_gs
from ModelPredictiveControl.MyMath import solve_lin_gs_with_Y
from ModelPredictiveControl.MyMath import solve_lin_gs_structured


class QuadraticProgram:

    def __init__(self, T, n, m):

        self.kappa = 90  # >0 barrier parameter

        self.T = T
        self.n = n
        self.m = m
        self.A = np.zeros([n, n])
        self.B = np.zeros([n, m])
        self.C = np.zeros([T*n, T*(n+m)])
        self.b = np.zeros([T*n, 1])

        self.H = np.eye(T*(n+m), T*(n+m))
        self.g = np.zeros([T*(m+n), 1])
        self.r = np.zeros([m, 1])
        self.S = np.zeros_like(self.B)

        self.Fx = None
        self.f = None
        # TODO Instance attributes als None initieren und abfragen, ob gesetzt

    def set_sys_dynamics(self, A, B):

        T, n, m = self.T, self.n, self.m

        self.A[:] = A
        self.B[:] = B

        # Equality constraints
        self.C[0:n, 0:m+n] = np.hstack([-B, np.eye(n, n)])
        for i in range(1, T):
            self.C[i*n:(i+1)*n, m+(i-1)*(m+n):m+i*(m+n)+n] = np.hstack([-A, -B, np.eye(n, n)])

        self.b = np.zeros([T*n, 1])  # TODO add disturbance if not zero

    def set_weighting(self, Q, q, R, r, S, Qf, qf):

        T, n, m = self.T, self.n, self.m

        self.r = r
        self.S = S

        # Cost function
        H = np.eye(T*(n+m), T*(n+m))
        H[0:m, 0:m] = R
        QSR = np.hstack([np.vstack([Q, S.T]), np.vstack([S, R])])
        for i in range(1, T):
            H[m+(i-1)*(m+n):m+i*(m+n), m+(i-1)*(m+n):m+i*(m+n)] = QSR
        H[m+(T-1)*(m+n):m+(T-1)*(m+n)+n, m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = Qf
        self.H = H

        g = np.zeros([T*(m+n), 1])
        for i in range(0, T):
            g[i*(n+m):(i+1)*(n+m)] = np.vstack([r, q])
        g[(T-1)*(n+m)+m:T*(n+m)] = qf
        self.g = g
        pass

    def set_constraints(self, Fu, fu, Fx, fx, Ff, ff):

        T, n, m = self.T, self.n, self.m

        self.Fx = Fx

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

        self.f = f

        h = np.zeros([T*np.shape(f)[0]+np.shape(ff)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f
        h[T*np.shape(f)[0]:T*np.shape(f)[0]+np.shape(ff)[0]] = ff
        self.h = h

    def h_of_xk(self, xk):
        h = np.zeros_like(self.h) + self.h  # does not change self.h
        h[0:np.shape(self.Fx)[0]] -= np.dot(self.Fx, xk)
        return h

    def g_of_xk(self, xk):
        pass

    def b_of_xk(self, xk):
        pass

    def form_d(self, xk, zv_k):
        # Form d for further use
        h = self.h_of_xk(xk)
        d = np.zeros([np.shape(self.P)[0], 1])
        d[:] = 1/(h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))
        return d

    def form_Phi(self, d):
        Phi = 2*self.H\
              + self.kappa*np.dot(np.dot(self.P.T, matrix_diag(d*d)), self.P)
        return Phi

    def solve(self, xk, zv_k):

        T = self.T
        n = self.n
        m = self.m

        self.b[0:self.n] = np.dot(self.A, xk)
        self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)

        d = self.form_d(xk, zv_k)
        Phi = self.form_Phi(d)

        r = np.vstack(self.residual(xk, zv_k))
        rd, rp = self.residual(xk, zv_k)
        SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.eye(self.C.shape[0], self.C.shape[0])*0])])

        v = zv_k[self.T*(self.m+self.n):]
        lsg = solve_lin_gs_structured(Phi, r, self.A, self.B, self.C, T, m, n, v)
        # lsg = solve_lin_gs_with_Y(Phi, self.C, rd, rp, v)
        # print(lsg)

        # lsg = solve_lin_gs(SS, -r)
        # print(lsg)
        return lsg

    def residual(self, xk, zv_k):

        d = self.form_d(xk, zv_k)
        rd = 2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):])
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - self.b

        if not self.check(xk, zv_k):
            return rd + 100000000000000000000000000000000000000, rp + 100000000000000000000000000000000000000

        return rd, rp

    def old_residual(self, xk, zv_k):
        h = np.zeros_like(self.h) + self.h
        h[0:np.shape(self.Fx)[0]] = self.f - np.dot(self.Fx, xk)

        d = np.zeros([np.shape(self.P)[0], 1])
        d[:] = 1/(h[:]-np.dot(self.P[:], zv_k[0:self.T*(self.m+self.n)]))

        # print 'term12',np.square(2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, d)).sum()
        rd = 2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)]) + self.g + self.kappa*np.dot(self.P.T, d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):])
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - self.b

        if not self.check(xk, zv_k):
            return rd + 100000000000000000000000000000000000000, rp + 100000000000000000000000000000000000000

        return rd, rp

    def check(self, xk, zv_k):
        h = self.h_of_xk(xk)
        return ((np.dot(self.P, zv_k[0:self.T*(self.m+self.n)]) - h) < 0).all()



