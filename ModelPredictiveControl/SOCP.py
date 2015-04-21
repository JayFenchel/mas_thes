#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from ModelPredictiveControl.MyMath import matrix_diag
from ModelPredictiveControl.MyMath import solve_lin_gs_structured

class SOCP:

    def __init__(self, T, n, m):

        # Kappa
        self.kappa = 90  # >0 barrier parameter
        # Dimensions
        self.T = T
        self.n = n
        self.m = m
        # System declaration
        self.A = np.zeros([n, n])
        self.B = np.zeros([n, m])
        self.C = np.zeros([T*n, T*(n+m)])
        self.b = np.zeros([T*n, 1])
        # Weighting for cost function
        self.H = np.eye(T*(n+m), T*(n+m))
        self.g = np.zeros([T*(m+n), 1])
        self.r = np.zeros([m, 1])
        self.S = np.zeros_like(self.B)
        # Linear constraints
        self.Fx = None
        self.f = None
        # allgemeine constraints
        self.P = None
        self.h = None
        # Reference
        self.u_ref = None
        self.z_ref = None
    # self.z_ref = np.zeros([T*(m+n), 1])  # for test case
        self.ref_update = None
    # # TODO
   #  self.Ff_qc = None
   #  self.socc_A = None

    def set_sys_dynamics(self, A, B):

        T, n, m = self.T, self.n, self.m

        self.A[:] = A
        self.B[:] = B

        # Equality constraints
        self.C[0:n, 0:m+n] = np.hstack([-B, np.eye(n, n)])
        for i in range(1, T):
            self.C[i*n:(i+1)*n, m+(i-1)*(m+n):m+i*(m+n)+n] = np.hstack([-A, -B, np.eye(n, n)])

        self.b = np.zeros([T*n, 1])  # TODO add disturbance if not zero

    def set_ref_trajectory(self, x_ref):
        # TODO Check dimensions
        # TODO Reference as a function
        # TODO Reference for Input u
        if self.z_ref is not None:
            self.z_ref = np.zeros([self.T*(self.m+self.n), 1])
        end = min(np.shape(x_ref)[1], self.T)  # Länge der vorgegebenen Reference
        self.u_ref = np.zeros([self.m, 1])
        for i in range(end):
            self.z_ref[i*(self.n+self.m):(i+1)*(self.n+self.m)] = np.vstack([self.u_ref, x_ref[:, i:i+1]])
        for i in range(end, self.T):  # wird nur ausgeführt, wenn T>Anzahl der Referenzpunkte
            self.z_ref[i*(self.n+self.m):(i+1)*(self.n+self.m)] = np.vstack([self.u_ref, x_ref[:, -1:]])
        # übrige Werte hinterlegen
        self.ref_update = x_ref[end:np.shape(x_ref)[1]]

    def update_ref_trajectory(self):
        self.z_ref[:-(self.n+self.m)] = self.z_ref[self.n+self.m:]
        if np.shape(self.ref_update)[0] > 0:
            self.z_ref[-(self.n+self.m):] = np.vstack([self.u_ref, self.ref_update[0:1].T])
            self.ref_update = np.delete(self.ref_update, 0, 0)

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

    def set_lin_constraints(self, Fu, fu, Fx, fx, Ff, ff):

        T, n, m = self.T, self.n, self.m

        self.Fx = Fx
        self.Fu = Fu
        self.Ff = Ff

        f = np.vstack([fx, fu]) # stacken - anderer branch
        f = fx + fu  # TODO fraglich

        self.f = f

        h = np.zeros([T*np.shape(f)[0]+np.shape(ff)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f
        h[T*np.shape(f)[0]:T*np.shape(f)[0]+np.shape(ff)[0]] = ff
        self.h = h

    # def add_qc(self, Ff_qc=None, alpha=None):
    #     # add quadratic constraint
    #     # TODO auf ausführliche Form, siehe Zettel erweitern
    #     self.Ff_qc = Ff_qc
    #     if alpha is not None:
    #         self.h = np.vstack([self.h, alpha])
    #
    # def add_socc(self, socc_A=None, socc_b=None, socc_c=None, socc_d=None):
    #     # add second-order cone constraint
    #     if socc_d is not None:
    #         self.h = np.vstack([self.h, socc_d*socc_d])
    #
    #     self.socc_A = socc_A
    #     self.socc_b = socc_b
    #     self.socc_c = socc_c
    #     self.socc_d = socc_d
    #
    #     pass
    #
    # def h_of_xk(self, xk):
    #     h = np.zeros_like(self.h) + self.h  # does not change self.h
    #     h[0:np.shape(self.Fx)[0]] -= np.dot(self.Fx, xk)
    #     return h
    #
    # def g_of_xk(self, xk):
    #     pass
    #
    # def b_of_xk(self, xk):
    #     pass
    #
    # def _A_of_socc_A_b(self, zk):
    #     # TODO [0:5] nur als Behelf, um die Dimensionen richtig zu machen
    #     _A_ = np.zeros([(self.m + self.n) * self.T, 1]) # Muss eine ganze Zeile für P werden
    #     _A_[-5:] = 2*(-self.socc_d*self.socc_c - np.dot(self.socc_c.T, zk[-5:])*self.socc_c + np.dot(self.socc_A.T, np.dot(self.socc_A, zk[-5:]) + self.socc_b))
    #     return _A_.T
    #
    # def d_A_dz_of_socc_A_b(self, zk):
    #     return 2*(-self.socc_c.T*self.socc_c.T + np.dot(self.socc_A.T, self.socc_A))

    def P_of_zk(self, zk):

        T, n, m = self.T, self.n, self.m
        # Inequality constraints
        n_Fu = np.shape(self.Fu)[0]
        P = np.zeros([T*n_Fu+np.shape(self.Ff)[0], T*(n+m)])
        P[0:n_Fu, 0:m] = self.Fu
        for i in range(1, T):
            Hilf = np.hstack([self.Fx, self.Fu])
            P[i*n_Fu:(i+1)*n_Fu, m+(i-1)*(m+n):m+i*(m+n)] = Hilf

        P[T*n_Fu:T*n_Fu+np.shape(self.Ff)[0], m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = self.Ff
        self.P = P
        P = np.zeros([np.shape(self.P)[0], np.shape(self.P)[1]])
        P += self.P  # does not change self.P
        if self.Ff_qc is not None:
            P = np.vstack([P, np.dot(zk.T, self.Ff_qc)[:]]) #  bei socc nur zk[-5:] genommen, um auf die Zeile in P zu kommen
        if self.socc_A is not None and self.socc_b is not None and self.socc_c is not None:
            _A_ = self._A_of_socc_A_b(zk)
            P = np.vstack([P, _A_])

        return P

    def form_d(self, xk, zv_k):
        # Form d for further use
        P = self.P_of_zk(zv_k[0:self.T*(self.m+self.n)])
        h = self.h_of_xk(xk)
        d = np.zeros([np.shape(P)[0], np.shape(zv_k)[1]])
        d[:] = 1/(h[:]-np.dot(P[:], zv_k[0:self.T*(self.m+self.n)]))
        return d

    def form_Phi(self, d, zk):
        P = self.P_of_zk(2*zk)
        if self.Ff_qc is not None:
            term_for_qc = d[-1]*2*self.Ff_qc
        else:
            term_for_qc = 0
        Phi = 2*self.H\
              + self.kappa*(np.dot(np.dot(P.T, matrix_diag(d*d)), P) + term_for_qc)  # *2 siehe Zettel
        return Phi

    def solve(self, xk, zv_k):

        T = self.T
        n = self.n
        m = self.m

        self.b[0:self.n] = np.dot(self.A, xk)
        self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)

        d = self.form_d(xk, zv_k)
        Phi = self.form_Phi(d, zv_k[0:self.T*(self.m+self.n)])

        rd, rp = self.residual(xk, zv_k)

        lsg = solve_lin_gs_structured(Phi, rd, rp, self.A, self.B, self.C, T, n, m)
        return lsg

    def residual_norm(self, zv_k, xk):

        if not self.check(xk, zv_k):
            return np.nan
        return np.square(np.vstack(self.residual(xk, zv_k))).sum(0)

    def residual(self, xk, zv_k):

        z_soll = np.zeros_like(zv_k[0:self.T*(self.m+self.n)])
        for i in range(self.m + 3, self.T*(self.m+self.n), self.m+self.n):
            z_soll[i] = 100
        # print(z_soll)
        d = self.form_d(xk, zv_k)
        rd = (2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)] - self.z_ref) + self.g
             + self.kappa*np.dot(self.P_of_zk(2*zv_k[0:self.T*(self.m+self.n)]).T, d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):]))
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - self.b
        return rd, rp

    def check(self, xk, zv_k):
        h = self.h_of_xk(xk)
        return ((np.dot(self.P_of_zk(zv_k[0:self.T*(self.m+self.n)]), zv_k[0:self.T*(self.m+self.n)]) - h) < 0).all()