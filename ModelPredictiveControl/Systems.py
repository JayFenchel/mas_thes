#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from numpy import diag


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
        self.delta_t = 1  # Länge der Schritte # TODO richtige Zeitschitte einbauen

class AirCraft:
    def __init__(self):

        self.T = 10
        self.delta_t = 1  # Länge der Schritte # TODO richtige Zeitschitte einbauen
        # mu = 100

        # discrete-time system
        Ad = np.array([[  0.23996015,   0., 0.17871287,   0., 0.],
                       [ -0.37221757,   1., 0.27026411,   0., 0.],
                       [ -0.99008755,   0., 0.13885973,   0., 0.],
                       [-48.93540655, 64.1, 2.39923411,   1., 0.],
                       [0., 0., 0., 0., 0.]])
        self.A = Ad
        Bd = np.array([[-1.2346445 ],
                       [-1.43828223],
                       [-4.48282454],
                       [-1.79989043],
                       [1.]])
        self.B = Bd
        n = Ad.shape[1]  # columns in A
        self.n = n
        m = Bd.shape[1]  # columns in B
        self.m = m

        # Weighting matrices for a problem with a better condition number
        self.Q = np.array(diag([1014.7, 3.2407, 5674.8, 0.3695, 471.75]))
        self.q = np.zeros([n, 1])
        self.R = np.array(diag([471.65]))
        self.r = np.zeros([m, 1])
        P = self.Q  # TODO Was ist P, terminal weighting?
        self.Qf = P
        qf = np.zeros([n, 1])
        self.qf = qf
        self.S = np.zeros_like(Bd)  # no combined weighting

        # input constraints
        eui = 0.262  # rad (15 degrees). Elevator angle.
        u_lb = -eui
        u_ub = eui
        Ku = np.array([[0.],
              [0.],
              [-1],
              [0],
              [0.],
              [0.],
              [1],
              [0],])
        self.Fu = Ku
        fu = np.zeros([np.shape(Ku)[0], 1])
        fu[np.shape(Ku)[0]/2-2] = -(u_lb)
        fu[np.shape(Ku)[0]/2-1] = -(u_lb)  # TODO wirklich mit den 0en?
        fu[np.shape(Ku)[0]-2] = (u_ub)
        fu[np.shape(Ku)[0]-1] = (u_ub)
        print fu
        # mixed constraints
        ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
        ex5 = 0.524 * self.delta_t  # rad/s * dt input slew rate constraint in discrete time
        ey3 = 30.
        # bounds
        e_lb = [[-ex2], [-ey3], [-ex5], [0]]
        e_ub = [[ex2], [ey3], [ex5], [0]]
        # constraint matrices
        Kx = np.array([[0, -1., 0., 0., 0.],
              [128.2, -128.2, 0, 0., 0.],
              [0., 0., 0., 0., 1.],
              [0., 0., 0., 0., 0.],
              [0., 1., 0., 0., 0.],
              [-128.2, 128.2, 0, 0., 0.],
              [0., 0., 0., 0., -1.],
              [0., 0., 0., 0., 0.]])
        self.Fx = Kx
        fx = np.ones([2*4, 1])
        fx [0:4] = - np.array(e_lb)
        fx[4:2*4] = np.array(e_ub)

        f = np.vstack([fx, fu]) # stacken - anderer branch
        f = fx + fu  # TODO fraglich
        self.f = f

        # terminal state constraints
        self.ff = fx
        self.ff[3] = 1
        self.ff[7] = 1
        self.Ff = Kx