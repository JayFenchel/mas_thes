#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np


class SimpleExample:
    def __init__(self):
        # System Dynamics and Control
        # TODO Disturbance w(t) einführen
        A = np.array([[-1, 2, 0],
                      [0, -2, 1],
                      [0, 1, 1]])
        self.A = A
        B = np.array([[0, 0],
                      [0, 0],
                      [1, 0]])
        self.B = B

        n = A.shape[1]  # columns in A
        self.n = n
        m = B.shape[1]  # columns in B
        self.m = m

        # Objective
        self.Q = np.eye(n, n)  # weighting of states
        self.q = np.eye(n, 1)*0
        self.R = np.eye(m, m)  # weighting of inputs
        self.r = np.eye(m, 1)*0
        self.S = np.zeros_like(B)  # non zero if not separable in states and control

        # Box constraints
        Fx = np.eye(2*(n+m), n)
        Fx[n:2*n] = Fx[0:n]
        self.Fx = Fx
        fx = np.ones([2*n, 1])
        fx[n:2*n] = -1
        fx *= 10
        Fu = np.eye(2*(n+m), m)
        Fu[2*n:2*n+m],Fu[2*n+m:2*(n+m)], Fu[0:m] = Fu[0:m], Fu[0:m], 0
        self.Fu = Fu
        fu = np.ones([2*m, 1])
        fu[m:2*m] = -1
        fu *= 2
        f = np.vstack([fx, fu])
        self.f = f

        # Terminal Cost, terminal state constraint
        # erstmall alles Null
        self.Qf = np.zeros_like(A)
        qf = np.eye(n, 1)
        qf[0] = 0
        self.qf = qf
        self.Ff = np.zeros_like(A)
        self.ff = np.zeros_like(qf)+1 # TODO, was wenn keine ff12

        self.T = 10  # Planning Horizon