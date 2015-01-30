#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np



class QuadraticProgram:

    def __init__(self, sys):

        n = sys.n
        m = sys.m
        T = sys.T  # Planning horizon

        self.n = n
        self.m = m
        self.T = T
        self.A = sys.A

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

        C = np.eye((T)*n, (T)*(n+m))
        C[0:n].T[0:m+n] = np.vstack([-sys.B.T, np.eye(n, n)])
        for i in range(1, T):
            C[i*n:(i+1)*n].T[m+(i-1)*(m+n):m+i*(m+n)+n] = np.vstack([-sys.A.T, -sys.B.T, np.eye(n, n)])
        self.C = C

        self.b = np.eye(T*n, 1)

    def solve(self, xk):

        T = self.T
        n = self.n
        m = self.m

        self.b[0:n] = np.dot(self.A, xk)




        return xk
