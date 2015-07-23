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
        self.roh = 10  # Wert aus paper
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
        # soft constraints
        self.P_soft = None
        self.h_soft = None
        self.Fx_soft = None
        self.f_soft = None
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

    def set_problem(self):
        T = self.T

        if self.h is None:
            print('No linear constraints have been set')
            self.h = np.zeros([0, 1])

        # add quadratic constraint to h
        if self.qc is not None:
            for qc in self.qc:
                h_add = np.zeros([T-1, 1])
                for i in range(0, T-1):
                    h_add[i] = qc['alpha']
                self.h = np.vstack([self.h, h_add])

        # add quadratic constraint (end) to h
        if self.qc_end is not None:
            for qc in self.qc_end:
                self.h = np.vstack([self.h, qc['alpha']])

        # add second-order cone constraint to h
        if self.socc is not None:
            for socc in self.socc:
                h_add = np.zeros([T-1, 1])
                for i in range(0, T-1):
                    h_add[i] = socc['d']*socc['d'] -\
                               np.dot(socc['b'].T, socc['b'])
                self.h = np.vstack([self.h, h_add])

        # add second-order cone constraint (end) to h
        if self.socc_end is not None:
            for socc in self.socc_end:
                self.h = np.vstack([self.h, socc['d']*socc['d'] -
                                    np.dot(socc['b'].T, socc['b'])])

        print('problem is set!')

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

        # h vector aus T-1 mal f und ff
        f = fx + fu

        h = np.zeros([T*np.shape(f)[0]+np.shape(ff)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f
        h[T*np.shape(f)[0]:T*np.shape(f)[0]+np.shape(ff)[0]] = ff
        self.h = h

        self.Fx = Fx
        self.f = f
        self.ff = ff

    def set_soft_constraints(self, Fu, fu, Fx, fx):
        # Set soft constraints with use of Kreisselmeier-Steinhauser Function

        T, n, m = self.T, self.n, self.m
        n_Fu = np.shape(Fu)[0]

        # P_tilde matrix with FuFx blocks
        P = np.zeros([T*n_Fu, T*(n+m)])
        P[0:n_Fu, 0:m] = Fu
        for i in range(1, T):
            Hilf = np.hstack([Fx, Fu])
            P[i*n_Fu:(i+1)*n_Fu, m+(i-1)*(m+n):m+i*(m+n)] = Hilf
        self.P_soft = P

        self.Fx_soft = Fx
        self.f_soft = fx + fu

    # Adding a quadratic constraint (type='end' for final constraints)
    def add_qc(self, gamma=None, beta=None, alpha=None, type='trajectory'):
        if beta is None:
            beta = np.zeros([np.shape(gamma)[0], 1])
        if gamma is None or alpha is None:
            print('gamma and alpha have to be not None.')

        if type == 'trajectory':
            if self.qc is not None:
                self.qc.append({'gamma': gamma, 'beta': beta, 'alpha': alpha})
            else:
                self.qc = [{'gamma': gamma, 'beta': beta, 'alpha': alpha}]
            pass

        elif type == 'end':
            if self.qc_end is not None:
                self.qc_end.append({'gamma': gamma, 'beta': beta, 'alpha': alpha})
            else:
                self.qc_end = [{'gamma': gamma, 'beta': beta, 'alpha': alpha}]

    # Adding a second order cone constraint (type='end' for final constraints)
    def add_socc(self, socc_A=None, socc_b=None,
                 socc_c=None, socc_d=None, type='trajectory'):
        if type == 'trajectory':
            if self.socc is not None:
                self.socc.append({'A': socc_A, 'b': socc_b, 'c': socc_c, 'd': socc_d})
            else:
                self.socc = [{'A': socc_A, 'b': socc_b, 'c': socc_c, 'd': socc_d}]

        elif type == 'end':
            if self.socc_end is not None:
                self.socc_end.append({'A': socc_A, 'b': socc_b, 'c': socc_c, 'd': socc_d})
            else:
                self.socc_end = [{'A': socc_A, 'b': socc_b, 'c': socc_c, 'd': socc_d}]
        else:
            print('SOCC: type is not known')
            exit()

    def calculate_kappa(self, zk):
        dim = self.n+self.m
        cost = np.dot(zk.T, np.dot(self.H, zk)) + np.dot(self.g.T, zk)
        kappa = 0.01*cost/dim
        return kappa

    def h_of_xk(self, xk):

        T, n, m = self.T, self.n, self.m

        if self.Fx is not None:
            h = np.zeros_like(self.h)
            h[:] = self.h[:]
            h[0:np.shape(self.Fx)[0]] -= np.dot(self.Fx, xk)
        else:
            h = np.zeros_like(self.h)
            h[:] = self.h[:]

        return h

    def h_soft_of_xk(self, xk):
        # for soft constraints
        T, n, m = self.T, self.n, self.m
        f = self.f_soft

        h = np.zeros([T*np.shape(f)[0], 1])
        for i in range(0, T):
            h[i*np.shape(f)[0]:(i+1)*np.shape(f)[0]] = f

        h[0:np.shape(self.Fx_soft)[0]] -= np.dot(self.Fx_soft, xk)

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
        # TODO Konstante Terme nur einmal berechnen
        p_i = np.dot(socc['A'].T, (np.dot(socc['A'], xk_) + 2*socc['b'])) -\
            socc['c']*(np.dot(socc['c'].T, xk_) + 2*socc['d'])
        return p_i.T

    def P_of_zk(self, zk):

        T, n, m = self.T, self.n, self.m
        # Inequality constraints
        if self.P is None:
            P = np.zeros([0, T*(n+m)])
        else:
            P = self.P

        # add quadratic constraint to P
        if self.qc is not None:
            for qc in self.qc:
                P_add = np.zeros([T-1, np.shape(P)[1]])
                for i in range(0, T-1):  # nur bis T-1, da T Index für qc_end
                    P_add[i, m+i*(m+n):m+i*(m+n)+n] =\
                        qc['beta'].T + np.dot(zk[m+i*(n+m):m+i*(n+m)+n].T, qc['gamma'])
                        # TODO obige Zeile als Funktion auslagern
                P = np.vstack([P, P_add])

        # add quadratic constraint (end) to P
        if self.qc_end is not None:
            for qc in self.qc_end:
                P_add = np.zeros([1, np.shape(P)[1]])
                P_add[0, m+(T-1)*(n+m):m+(T-1)*(n+m)+n] =\
                    qc['beta'].T + np.dot(zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n].T, qc['gamma'])
                    # TODO obige Zeile als Funktion auslagern
                P = np.vstack([P, P_add])

        # add second-order cone constraint to P
        if self.socc is not None:
            for socc in self.socc:
                P_add = np.zeros([T-1, np.shape(P)[1]])
                for i in range(0, T-1):  # nur bis T-1, da T Index für socc_end
                    P_add[i, m+i*(m+n):m+i*(m+n)+n] =\
                        self._A_of_socc(socc, zk[m+i*(n+m):m+i*(n+m)+n])
                P = np.vstack([P, P_add])

        # add second-order cone constraint (end) to P
        if self.socc_end is not None:
            for socc in self.socc_end:
                P_add = np.zeros([1, np.shape(P)[1]])
                P_add[0, m+(T-1)*(n+m):m+(T-1)*(n+m)+n] =\
                    self._A_of_socc(socc, zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n])
                P = np.vstack([P, P_add])

        return P

    def form_d(self, xk, zv_k):
        T, n, m = self.T, self.n, self.m
        # Form d for further use
        P = self.P_of_zk(zv_k[0:self.T*(self.m+self.n)])
        h = self.h_of_xk(xk)
        d = np.zeros([np.shape(P)[0], np.shape(zv_k)[1]])
        d[:] = 1/(h[:]-np.dot(P[:], zv_k[0:self.T*(self.m+self.n)]))

        return d

    def form_d_soft(self, xk, zv_k):
        P = self.P_soft
        h = self.h_soft_of_xk(xk)
        d = np.zeros([np.shape(P)[0], np.shape(zv_k)[1]])
        e_term = np.exp(self.roh*(h[:]-np.dot(P[:], zv_k[0:self.T*(self.m+self.n)])))
        d[:] = 1/(1 + e_term)
        return d

    def form_d_soft_dach(self, xk, zv_k):
        P = self.P_soft
        h = self.h_soft_of_xk(xk)
        d = np.zeros([np.shape(P)[0], np.shape(zv_k)[1]])
        e_term = np.exp(self.roh*(h[:]-np.dot(P[:], zv_k[0:self.T*(self.m+self.n)])))
        d[:] = e_term*1/((1 + e_term)*(1 + e_term))
        return d

    def term_for_qc(self, zk):
        T, n, m = self.T, self.n, self.m
        # add term for qc
        term_for_qc = np.zeros([T*(n+m), T*(n+m)])
        if self.qc is not None:
            for qc in self.qc:
                for i in range(0, T-1):  # nur bis T-1, da T Index für qc_end
                    d_k = 1/(qc['alpha'] - np.dot(
                        qc['beta'].T + np.dot(zk[m+i*(n+m):m+i*(n+m)+n].T, qc['gamma']),
                        zk[m+i*(n+m):m+i*(n+m)+n]))
                    # nur passender  nxn-Block für jeweile (T-1) x_k
                    term_for_qc[m+i*(n+m):m+i*(n+m)+n, m+i*(n+m):m+i*(n+m)+n] +=\
                        d_k*2*qc['gamma']
        # add term for qc (end)
        term_for_qc_end = np.zeros([T*(n+m), T*(n+m)])
        if self.qc_end is not None:
            for qc in self.qc_end:
                d_k = 1/(qc['alpha'] - np.dot(
                    qc['beta'].T + np.dot(zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n].T, qc['gamma']),
                    zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n]))
                term_for_qc_end[m+(T-1)*(n+m):m+(T-1)*(n+m)+n, m+(T-1)*(n+m):m+(T-1)*(n+m)+n] +=\
                    d_k*2*qc['gamma']  # nur unterer rechter Block nxn bei end_qc

        return term_for_qc + term_for_qc_end

    def term_for_socc(self, zk):
        T, n, m = self.T, self.n, self.m
        # add term for qc
        term_for_socc = np.zeros([T*(n+m), T*(n+m)])
        if self.socc is not None:
            for socc in self.socc:
                for i in range(0, T-1):  # nur bis T-1, da T Index für qc_end
                    d_k = 1/((np.dot(socc['c'].T, zk[m+i*(n+m):m+i*(n+m)+n]) + socc['d'])*(np.dot(socc['c'].T, zk[m+i*(n+m):m+i*(n+m)+n]) + socc['d']) -
                             ((np.dot(socc['A'], zk[m+i*(n+m):m+i*(n+m)+n])+socc['b'])* (np.dot(socc['A'], zk[m+i*(n+m):m+i*(n+m)+n])+socc['b'])).sum())
                    # nur passender  nxn-Block für jeweilige (T-1) x_k
                    term_for_socc[m+i*(n+m):m+i*(n+m)+n, m+i*(n+m):m+i*(n+m)+n] +=\
                        d_k*-2*(np.dot(socc['c'], socc['c'].T) - np.dot(socc['A'].T, socc['A']))
        # add term for qc (end)
        term_for_socc_end = np.zeros([T*(n+m), T*(n+m)])
        if self.socc_end is not None:
            for socc in self.socc_end:
                d_k = 1/((np.dot(socc['c'].T, zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n]) + socc['d'])*(np.dot(socc['c'].T, zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n]) + socc['d']) -
                         ((np.dot(socc['A'], zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n])+socc['b'])*(np.dot(socc['A'], zk[m+(T-1)*(n+m):m+(T-1)*(n+m)+n])+socc['b'])).sum())
                term_for_socc_end[m+(T-1)*(n+m):m+(T-1)*(n+m)+n, m+(T-1)*(n+m):m+(T-1)*(n+m)+n] +=\
                    d_k*-2*(np.dot(socc['c'], socc['c'].T) - np.dot(socc['A'].T, socc['A']))  # nur unterer rechter Block nxn bei end_qc

        return term_for_socc + term_for_socc_end

    def form_Phi(self, d, zk):

        P = self.P_of_zk(2*zk)
        #TODO add Term for socp
        Phi = 2*self.H\
              + self.kappa*(np.dot(np.dot(P.T, matrix_diag(d*d)), P) +
                            self.term_for_qc(zk) + self.term_for_socc(zk))  # *2 siehe Zettel
        return Phi

    def form_Phi_soft(self, d, zk):
        Phi = self.roh*(np.dot(np.dot(self.P_soft.T, matrix_diag(d)), self.P_soft))
        return Phi

    def solve(self, xk, zv_k):

        T = self.T
        n = self.n
        m = self.m

        self.g[0:m] = self.r + 2*np.dot(self.S.T, xk)

        d = self.form_d(xk, zv_k)
        Phi = self.form_Phi(d, zv_k[0:self.T*(self.m+self.n)])
        if self.P_soft is not None:
            d_soft_dach = self.form_d_soft_dach(xk, zv_k)
            Phi += self.form_Phi_soft(d_soft_dach, zv_k[0:self.T*(self.m+self.n)])

        rd, rp = self.residual(xk, zv_k)

        # SS = np.hstack([np.vstack([Phi, self.C]), np.vstack([self.C.T, np.zeros([self.C.shape[0], self.C.shape[0]])])])
        # lsg = linalg.solve(SS, -np.vstack([rd, rp]))

        # q, r = householder(SS) # TODO housholder trafo scheint hier nicht richtig zu funktionieren -> Test schreiben
        # TODO Ausgabe bei Div durch 0 in housholder
        # lsg1 = backward_substitution(r, np.dot(q.T, -np.vstack([rd, rp])))
        lsg = solve_lin_gs_structured(Phi, rd, rp, self.A, self.B, self.C, T, n, m, reg=0.00000001)
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
        if self.P_soft is not None:  # if there are defined soft constraints
            d_soft = self.form_d_soft(xk, zv_k)
            rd += np.dot(self.P_soft.T, d_soft)
        return rd, rp

    def check(self, xk, zv_k):
        h = self.h_of_xk(xk)
        return ((np.dot(self.P_of_zk(zv_k[0:self.T*(self.m+self.n)]), zv_k[0:self.T*(self.m+self.n)]) - h) < 0).all()