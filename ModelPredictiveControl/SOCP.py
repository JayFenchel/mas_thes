#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
# from scipy import linalg
from ModelPredictiveControl.MyMath import matrix_diag
from ModelPredictiveControl.MyMath import solve_lin_gs_structured
from ModelPredictiveControl.MyMath import householder
from ModelPredictiveControl.MyMath import backward_substitution

class SOCP:

    def __init__(self, T, n, m, x0=None, u0=None):

        # Kappa
        self.kappa = 90  # >0 barrier parameter
        # Dimensions
        self.T = T
        self.n = n
        self.m = m
        self.x0 = x0
        self.u0 = u0
        # System declaration
        self.A = np.zeros([n, n])
        self.B = np.zeros([n, m])
        self.C = np.zeros([T*n, T*(n+m)])
        self.b = None
        # Weighting for cost function
        self.H = np.eye(T*(n+m), T*(n+m))
        self.g = np.zeros([T*(m+n), 1])
        self.r = np.zeros([m, 1])
        self.S = np.zeros_like(self.B)
        # Linear constraints
        self.Fx = None
        self.f = None
        self.ff = None
        # allgemeine constraints
        self.P = None
        self.h = None
        # Reference
        self.u_ref = None
        self.z_ref = None
        self.z_ref = np.zeros([T*(m+n), 1])  # for test case
        self.ref_update = None
        # Nonlinear constraints
        self.qc = None
        self.qc_end = None
        self.socc_end = None
        self.socc = None

    def set_sys_dynamics(self, A, B):

        T, n, m = self.T, self.n, self.m

        self.A[:] = A
        self.B[:] = B

        # Equality constraints
        self.C[0:n, 0:m+n] = np.hstack([-B, np.eye(n, n)])
        for i in range(1, T):
            self.C[i*n:(i+1)*n, m+(i-1)*(m+n):m+i*(m+n)+n] = np.hstack([-A, -B, np.eye(n, n)])

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

    def set_weighting(self, Q=None, q=None, R=None, r=None, S=None,
                      Qf=None, qf=None):

        T, n, m = self.T, self.n, self.m
        # if None: set_zero(right_dimension)
        if q is None:
            q = np.zeros([n, 1])
        if r is None:
            r = np.zeros([m, 1])
        if S is None:
            S = np.zeros([n, m])
        if Qf is None and qf is None:
            Qf, qf = Q, q
        elif Qf is None or qf is None:
            if (q == np.zeros([n, 1])).all():
                qf = q
            else:
                print('Qf XOR qf is None, determinate qf = q or qf = 0')
        if Q is None or R is None or Qf is None:
            print('Some important weighting matrices are not defined!')
            exit()

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
        n_Fu = np.shape(Fu)[0]

        # P matrix with FuFx blocks
        P = np.zeros([T*n_Fu+np.shape(Ff)[0], T*(n+m)])
        P[0:n_Fu, 0:m] = Fu
        for i in range(1, T):
            Hilf = np.hstack([Fx, Fu])
            P[i*n_Fu:(i+1)*n_Fu, m+(i-1)*(m+n):m+i*(m+n)] = Hilf
        P[T*n_Fu:T*n_Fu+np.shape(Ff)[0], m+(T-1)*(m+n):m+(T-1)*(m+n)+n] = Ff
        self.P = P

        self.Fx = Fx

        self.f = fx + fu
        self.ff = ff

    # Adding a quadratic constraint (type='end' for final constraints)
    def add_qc(self, type='trajectory', F_qc=None, alpha=None):
        # TODO auf ausführliche Form, siehe Zettel erweitern
        if type == 'trajectory':
            if self.qc is not None:
                self.qc.append([F_qc, alpha])
            else:
                self.qc = [[F_qc, alpha]]
            pass

        elif type =='end':
            if self.qc_end is not None:
                self.qc_end.append([F_qc, alpha])
            else:
                self.qc_end = [[F_qc, alpha]]

    # Adding a second order cone constraint (type='end' for final constraints)
    def add_socc(self, type='trajectory', socc_A=None, socc_b=None,
                 socc_c=None, socc_d=None):
        if type == 'trajectory':
            if self.socc is not None:
                self.socc.append([socc_A, socc_b, socc_c, socc_d])
            else:
                self.socc = [[socc_A, socc_b, socc_c, socc_d]]

        elif type == 'end':
            if self.socc_end is not None:
                self.socc_end.append([socc_A, socc_b, socc_c, socc_d])
            else:
                self.socc_end = [[socc_A, socc_b, socc_c, socc_d]]
        else:
            print('SOCC: type is not known')
            exit()

    def h_of_xk(self, xk):

        T, n, m = self.T, self.n, self.m
        f = self.f

        h = np.zeros([T*np.shape(f)[0]+np.shape(self.ff)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f
        h[T*np.shape(f)[0]:T*np.shape(f)[0]+np.shape(self.ff)[0]] = self.ff

        h[0:np.shape(self.Fx)[0]] -= np.dot(self.Fx, xk)
        # add quadratic constraint (end) to h
        if self.qc_end is not None:
            for qc in self.qc_end:
                h = np.vstack([h, qc[1]])

        # add second-order cone constraint to h
        if self.socc is not None:
            for socc in self.socc:
                h_add = np.zeros([T-1, 1])
                for i in range(0, T-1):
                    h_add[i] = socc[3]*socc[3]  # TODO d²?
                h = np.vstack([h, h_add])

        # add second-order cone constraint (end) to h
        if self.socc_end is not None:
            for socc in self.socc_end:
                h = np.vstack([h, socc[3]*socc[3]])
        return h

    def g_of_xk(self, xk):
        pass

    def b_of_xk(self, xk):

        T, n = self.T, self.n
        b = np.zeros([T*n, 1])  # TODO add disturbance if not zero
        b[0:n] = np.dot(self.A, xk)
        return b

    # new line in P matrix corresponding to a second order cone constraint
    def _A_of_socc(self, socc, xk_):
        # TODO test _A_of_socc
        # TODO Konstante Terme nur einmal berechnen
        socc_A, socc_b, socc_c, socc_d = socc[0], socc[1], socc[2], socc[3]
        # TODO [0:n] nur als Behelf, um die Dimensionen richtig zu machen eigentlich jede Zeile mit den richtigen Einträgen aus z multiplizieren
        _A_ = 2*(-socc_d*socc_c - np.dot(socc_c.T, xk_)*socc_c
                 + np.dot(socc_A.T, np.dot(socc_A, xk_) + socc_b))
        return _A_.T

    def _A_of_qc(self, qc, zk):
        # return 2*z.T*Ff_alpha
        # TODO richtige Berechnung, richtige angabe in system.py
        pass

    def d_A_dz_of_socc_A_b(self, zk):
        return 2*(-self.socc_c_end.T*self.socc_c_end.T + np.dot(self.socc_A_end.T, self.socc_A_end))

    def P_of_zk(self, zk):

        T, n, m = self.T, self.n, self.m
        # Inequality constraints
        P = self.P

        # add quadratic constraint (end) to P
        if self.qc_end is not None:
            for qc in self.qc_end:
                P_add = np.dot(zk.T, qc[0])
                P = np.vstack([P, P_add])

        # add second-order cone constraint to P
        if self.socc is not None:
            for socc in self.socc:
                P_add = np.zeros([T-1, np.shape(self.P)[1]])
                for i in range(0, T-1):  # nur bis T-1, da T Index für socc_end
                    P_add[i, m+i*(m+n):m+i*(m+n)+n] =\
                        self._A_of_socc(socc, zk[m+i*(n+m):m+i*(n+m)+n])
                P = np.vstack([P, P_add])

        # add second-order cone constraint (end) to P
        if self.socc_end is not None:
            for socc in self.socc_end:
                P_add = np.zeros([1, np.shape(self.P)[1]])
                P_add[0, m+(T-1)*(n+m):m+(T-1)*(n+m)+n] =\
                    self._A_of_socc(socc, zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n])
                P = np.vstack([P, P_add])

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
        #TODO add Term for socp
        # add term for qc (end)
        term_for_qc = 0
        if self.qc_end is not None:
            for qc in self.qc_end:
                term_for_qc += d[-1]*2*qc[0]  # TODO Wirklich Summe über bilden?

        Phi = 2*self.H\
              + self.kappa*(np.dot(np.dot(P.T, matrix_diag(d*d)), P) + term_for_qc)  # *2 siehe Zettel
        return Phi

    def solve(self, xk, zv_k):

        T = self.T
        n = self.n
        m = self.m

        self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)

        d = self.form_d(xk, zv_k)
        Phi = self.form_Phi(d, zv_k[0:self.T*(self.m+self.n)])

        rd, rp = self.residual(xk, zv_k)

        # SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.zeros([self.C.shape[0], self.C.shape[0]])])])
        # lsg = linalg.solve(SS, -np.vstack([rd, rp]))

        # q, r = householder(SS) # TODO housholder trafo scheint hier nicht richtig zu funktionieren -> Test schreiben
        # TODO Ausgabe bei Div durch 0 in housholder
        # lsg1 = backward_substitution(r, np.dot(q.T, -np.vstack([rd, rp])))
        lsg = solve_lin_gs_structured(Phi, rd, rp, self.A, self.B, self.C, T, n, m, reg=0.000001)
        return lsg

    def residual_norm(self, zv_k, xk):

        if not self.check(xk, zv_k):
            return np.nan
        return np.square(np.vstack(self.residual(xk, zv_k))).sum(0)

    def residual(self, xk, zv_k):

        b = self.b_of_xk(xk)
        d = self.form_d(xk, zv_k)
        rd = (2*np.dot(self.H, zv_k[0:self.T*(self.m+self.n)] - self.z_ref) + self.g
             + self.kappa*np.dot(self.P_of_zk(2*zv_k[0:self.T*(self.m+self.n)]).T, d) + np.dot(self.C.T, zv_k[self.T*(self.m+self.n):]))
        rp = np.dot(self.C, zv_k[0:self.T*(self.m+self.n)]) - b
        return rd, rp

    def check(self, xk, zv_k):
        h = self.h_of_xk(xk)
        return ((np.dot(self.P_of_zk(zv_k[0:self.T*(self.m+self.n)]), zv_k[0:self.T*(self.m+self.n)]) - h) < 1e-08).all()