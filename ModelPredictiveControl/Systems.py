#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np
from scipy import io
from scipy.linalg import pinv
from numpy import diag
from ModelPredictiveControl.SOCP import SOCP
from ModelPredictiveControl.getdata import fromfile

test_dir = 'data/QP-Test-Problems/MAT_Files/'

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


# class Systems:
    # def __init__(self):
    #     self.QP = QuadraticProgram

def reorder(a, T, n, m):
    b = np.zeros_like(a)
    for i in range(0, T):
        b[:, i*(n+m):i*(n+m)+m] = a[:, n*T+i*m:n*T+(i+1)*m]
        b[:, i*(n+m)+m:i*(n+m)+m+n] = a[:, i*n:(i+1)*n]
    return b

def qp_from_test():
    epsilon = 1e-06  # um GNB als UNB zu schreiben
    # read data
    data = io.loadmat('%sCVXQP2_S.mat' % test_dir)
    # bounds
    low_bounds = data['lb']  # lower bounds
    up_bounds = data['ub']  # upper bounds
    # mixed inequality- and equality constraints of the form b_l <= A*x <= b_u
    b_l = data['rl']
    b_u = data['ru']
    A_hat = data['A'].toarray()

    n, m = np.shape(A_hat)
    T = 1

    if not (b_l == b_u).any():
        print('There are ONLY mixed inequality constraints')
        # kein epsilon verwenden
        exit()
    elif not (b_l == b_u).all():
        print('There are mixed inequality- and equality constraints')
        # epsilon verwenden
        b_l = b_l-epsilon
        b_u = b_u+epsilon
        exit()
    else:
        print('There are ONLY equality constraints')
        if not (b_l == b_u).all():
            print('Irgendetwas ist vorgefallen')
        b_hat = b_l
        x0 = b_hat
        # TODO Fall betrachten in denen nur einzelne Elemente Bounds = Inf
        u0 = np.zeros([m, 1])
        if (low_bounds > -np.inf).all():
            if (up_bounds < np.inf).all():
                u0[:] = .5*(low_bounds[:] + up_bounds[:])
            else:
                u0[:] = low_bounds[:] + 10
        else:
            if (up_bounds < np.inf).all():
                u0[:] = up_bounds[:] - 10
            else:
                u0 = np.zeros([m, 1])


        #transformed equality constraint of the form:
        #         x(t+1) = Ad*x(t) + Bd*u
        # m=0:    x(t+1) = Ad*x(t)
        # mit Ad = inv(A)*b*inv(x(t))

        Ad = np.zeros([n, n])
        Bd = A_hat
        qp = SOCP(T, n, m, x0=x0, u0=u0)
        qp.set_sys_dynamics(np.array(Ad), np.array(Bd))

        # weighting matrices
        Q = np.zeros([n, n])
        R = .5*data['Q'].toarray()
        r = data['c']
        qp.set_weighting(Q=Q, R=R, r=r)

        # bounds
        if not np.shape(low_bounds)[0] == m or not np.shape(up_bounds)[0] == m:
            print('Dimension of bounds are not equal to number of variables')
            exit()
        else:
            # input
            fu = np.zeros([2*(m+n), 1])
            fu[n:n+m] = -np.array(low_bounds)
            fu[2*n+m:2*(n+m)] = np.array(up_bounds)

            Ku = np.zeros([2*(m+n), m])
            Ku[n:n+m, :] = -np.eye(m)
            Ku[2*n+m:2*(n+m), :] = np.eye(m)

            # state
            fx = np.zeros([2*(m+n), 1])
            fx[0:n] = -(x0-0.00001)
            fx[n+m:2*n+m] = x0+0.00001

            Kx = np.zeros([2*(m+n), n])
            Kx[0:n, :] = -np.eye(n)
            Kx[n+m:2*n+m, :] = np.eye(n)

            ff = np.zeros([2*(n), 1])
            ff[0:n] = -(x0-0.00001)
            ff[n:2*n] = x0+0.00001

            Kf = np.zeros([2*(n), n])
            Kf[0:n, :] = -np.eye(n)
            Kf[n:2*n, :] = np.eye(n)

            qp.set_lin_constraints(Fu=Ku, fu=fu, Fx=Kx, fx=fx, Ff=Kf, ff=ff)


    return (qp, Bd, 2*x0)

def qp_from_new_sys():
    mm = io.loadmat('data/data_matrix.mat')  # load system matrices
    socp = io.loadmat('data/socp_matrices.mat')

    # discrete-time system
    Ad = mm['Asys']
    Bd = mm['Bsys']
    (n, m) = np.array(Bd).shape  # system dimensions
    T = 5  # prediction horizon
    dt = 0.5
    x0 = -socp['X0'][0:30]
    u0 = np.zeros([m, 1])
    qp = SOCP(T, n, m, x0=x0, u0=u0)
    qp.set_sys_dynamics(np.array(Ad), np.array(Bd))

    # weighting matrices
    Q = mm['Q_total']
    R = mm['R_total']
    P = Q  # terminal weighting
    qp.set_weighting(Q=Q, R=R, Qf=P)

    # input constraints
    eui = 0.262  # rad (15 degrees). Elevator angle.
    u_lb = -eui
    u_ub = eui
    Ku = np.array([[0.],
                   [0.],
                   [-1.],
                   [0.],
                   [0.],
                   [1.],])

    fu = np.zeros([np.shape(Ku)[0], 1])
    fu[np.shape(Ku)[0]/2-1] = -(u_lb)
    fu[np.shape(Ku)[0]-1] = (u_ub)

    # mixed constraints
    ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
    ex5 = 0.524 * dt  # rad/s * dt input slew rate constraint in discrete time
    ey3 = 30.

    # bounds
    e_lb = [[-ex2], [-ex5], [0]]
    e_ub = [[ex2], [ex5], [0]]
    fx = np.ones([2*3, 1])
    fx[0:3] = - np.array(e_lb)
    fx[3:2*3] = np.array(e_ub)

    (ncx, dummy) = np.array(e_ub).shape
    # constraint matrices
    Kx = np.zeros((ncx*2, n))
    # Kx[0, 6] = -1
    # Kx[3, 6] = 1
    # Kx[1, 24] = 1
    # Kx[4, 24] = -1

    # keine Platzhalterzeilen für u bounds in den end constraints
    Kf = np.zeros([4, n])
    ff = np.zeros([4, 1])
    # ff[:] = fx[:]
    ff[0:2] = fx[0:2]
    ff[2:4] = fx[3:5]  #TODO Werte raus aus SOCP file

    qp.set_lin_constraints(Fu=Ku, fu=fu, Fx=Kx, fx=fx, Ff=Kf, ff=ff)

    x_ref = np.array([[0.], [0.], [0.], [0.], [0.], [0.],
                      [0.], [0.], [0.], [0.], [0.], [0.],
                      [0.], [0.], [0.], [0.], [0.], [0.],
                      [200.], [0.], [0.], [0.], [0.], [0.],
                      [0.], [0.], [0.], [0.], [0.], [0.]])
    qp.set_ref_trajectory(x_ref)

    # V21 = reorder(socp['V21'][30:60, 30:], T, n, m)
    # V22 = reorder(socp['V22'][60:90, 30:], T, n, m)
    # V23 = reorder(socp['V22'][90:120, 30:], T, n, m)
    # V24 = reorder(socp['V22'][120:150, 30:], T, n, m)
    # V25 = reorder(socp['V22'][150:180, 30:], T, n, m)
    help = socp['V21']
    V21 = help[30:60, 30:60]

    help = socp['M21']
    M21 = np.array([help[0, 30:60]])
    c = socp['c_max']

    qp.add_socc(socc_A=np.array(V21), socc_c=M21.T,
                socc_b=np.zeros_like(M21.T), socc_d=c)
    qp.add_socc(type='end', socc_A=np.array(V21), socc_c=M21.T,
                socc_b=np.zeros_like(M21.T), socc_d=c)
    # TODO add same socc as type 'end'
    return qp

def qp_from_sys():
    # discrete-time system
    Ad = np.array([[  0.23996015,   0., 0.17871287,   0., 0.],
                   [ -0.37221757,   1., 0.27026411,   0., 0.],
                   [ -0.99008755,   0., 0.13885973,   0., 0.],
                   [-48.93540655, 64.1, 2.39923411,   1., 0.],
                   [0., 0., 0., 0., 0.]])
    A = Ad
    Bd = np.array([[-1.2346445 ],
                   [-1.43828223],
                   [-4.48282454],
                   [-1.79989043],
                   [1.]])
    B = Bd
    n = Ad.shape[1]  # columns in A
    m = Bd.shape[1]  # columns in B
    T = 10  # Prädiktionshorizont
    qp = SOCP(T, n, m)

    delta_t = 0.5

    qp.set_sys_dynamics(A, B)

    # Weighting matrices for a problem with a better condition number
    Q = np.array(diag([1014.7, 3.2407, 5674.8, 0.3695, 471.75]))
    q = np.zeros([n, 1])
    R = np.array(diag([471.65]))
    r = np.zeros([m, 1])
    P = Q  # terminal weighting
    qf = np.zeros([n, 1])
    S = np.zeros_like(Bd)  # no combined weighting

    qp.set_weighting(Q=Q, q=q, R=R, r=r, S=S, Qf=P, qf=qf)

    # input constraints
    eui = 0.262  # rad (15 degrees). Elevator angle.
    u_lb = -eui
    u_ub = eui
    Ku = np.array([[0.],
          [0.],
          [-1.],
          [0.],
          [0.],
          [0.],
          [1.],
          [0.],])
    Fu = Ku
    fu = np.zeros([np.shape(Ku)[0], 1])
    fu[np.shape(Ku)[0]/2-2] = -(u_lb)
    fu[np.shape(Ku)[0]/2-1] = -(u_lb)  # TODO wirklich mit den 0en?
    fu[np.shape(Ku)[0]-2] = (u_ub)
    fu[np.shape(Ku)[0]-1] = (u_ub)
    # mixed constraints
    ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
    ex5 = 0.524 * delta_t  # rad/s * dt input slew rate constraint in discrete time
    ey3 = 30.
    # bounds
    e_lb = [[-ex2], [-ey3], [-ex5], [0]]
    e_ub = [[ex2], [ey3], [ex5], [0]]
    # constraint matrices
    # TODO Matrizen ohne 0 Zeilen, um zu den Input constraints zu passen
    Kx = np.array([[0, -1., 0., 0., 0.],
          [128.2, -128.2, 0, 0., 0.],
          [0., 0., 0., 0., 1.],
          [0., 0., 0., 0., 0.],
          [0., 1., 0., 0., 0.],
          [-128.2, 128.2, 0, 0., 0.],
          [0., 0., 0., 0., -1.],
          [0., 0., 0., 0., 0.]])
    Fx = Kx
    fx = np.ones([2*4, 1])
    fx [0:4] = - np.array(e_lb)
    fx[4:2*4] = np.array(e_ub)

    # # terminal state constraints
    ff = fx
    ff[3] = 1
    ff[7] = 1
    Ff = Kx

    Ff_qc = np.zeros([T*(n+m), T*(n+m)])
    alpha = 1

    qp.set_lin_constraints(Fu, fu, Fx, fx, Ff, ff)
    qp.add_qc(type='end', F_qc=Ff_qc, alpha=alpha)
    x_ref = np.array([[0], [0], [0], [200], [0]])
    qp.set_ref_trajectory(x_ref)

    socc_A = np.array([[1, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0]])
    socc_b = np.array([[1],
                       [1]])
    socc_c = np.array([[1],
                      [0],
                      [0],
                      [0],
                      [0]])
    socc_d = np.array([[3]])

    qp.add_socc(socc_A=socc_A, socc_b=socc_b, socc_c=socc_c, socc_d=socc_d)
    qp.add_socc(socc_A=socc_A, socc_b=socc_b, socc_c=socc_c, socc_d=socc_d)
    qp.add_socc(type='end', socc_A=socc_A, socc_b=socc_b, socc_c=socc_c,
                socc_d=socc_d)
    qp.add_socc(type='end', socc_A=socc_A, socc_b=socc_b, socc_c=socc_c,
                socc_d=socc_d)
    return qp

class AirCraft:
    def __init__(self):

        self.T = 10
        self.delta_t = .5  # Länge der Schritte # TODO richtige Zeitschitte einbauen
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
              [-1.],
              [0.],
              [0.],
              [0.],
              [1.],
              [0.],])
        self.Fu = Ku
        fu = np.zeros([np.shape(Ku)[0], 1])
        fu[np.shape(Ku)[0]/2-2] = -(u_lb)
        fu[np.shape(Ku)[0]/2-1] = -(u_lb)  # TODO wirklich mit den 0en?
        fu[np.shape(Ku)[0]-2] = (u_ub)
        fu[np.shape(Ku)[0]-1] = (u_ub)
        print(fu)
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
