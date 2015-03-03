#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np


def vector_norm(a):
    n = 0
    for i in range(0, a.shape[0]):
        n += a[i]**2
    return np.sqrt(n)

def matrix_diag(a):

    return np.diag(a.T[0])

def householder(a):
    # Q-R-Zerlegung via Householdertransformation (nach Vorlage von Sager)
    neta, npe = a.shape  # neta - # of Zeilen, npe - # of Spalten in A
    Q = np.array(np.eye(neta))
    R = a

    for i in range(0, npe):
        w = np.array([R[i:neta, i]]).T
        alpha = vector_norm(w)
        v = (w - alpha*np.array(np.eye(neta-i, 1)))
        if vector_norm(v) != 0:
            v = v / np.linalg.norm(v)
        # else:
        #     print("Division durch 0 aufgetreten")

        H = np.array(np.eye(neta, neta))
        H[i:neta, i:neta] = np.array(np.eye(neta-i)) - 2*np.dot(v, v.T)
        Q = np.dot(H, Q)
        R = np.dot(H, R)

    return Q, R

def cholesky(hilf):
    a = np.zeros_like(hilf)
    a[:, :] = hilf  # damit die ursprüngliche Matrix nicht mitverändert wird
    # TODO Prüfen, ob a symmetrisch, A == A.T ?

    n = a.shape[0]
    for k in range(0, n):
        for j in range(0, k):
            a[k, k] = a[k, k] - a[k, j]*a[k, j]
        a[k, k] = np.sqrt(np.abs(a[k, k]))
        for i in range(k+1, n):
            for j in range(0, k):
                a[i, k] = a[i, k] - a[i, j]*a[k, j]
            a[i, k] = a[i, k]/a[k, k]

    for k in range(0, n):
        for j in range(k+1, n):
            a[k, j] = 0

    return(a)

def solve_lin_gs(A, b):

    neta, npe = A.shape  # neta - # of Zeilen, npe - # of Spalten in A
    Q, R = householder(A)
    # hhhh = cholesky(A)
    c = np.dot(Q, b)
    Rtilde = R[0:npe, 0:npe]
    ctilde = c[0:npe]
    x = backward_substitution(Rtilde, ctilde)
    return x

def form_Y(L_Phi, A, B, T, n, m):

    # TODO find ungenauigkeit, irgendwo ist vllt noch ein kleiner fehler
    Y = np.zeros([T*n, T*n])

    for i in range(0, T):
        for j in range(i, i+2):
            # Y[0, 0]
            if i == 0 and j == 0:
                Y[0:n, 0:n] = np.dot(B, backward_substitution(L_Phi.T[0:m, 0:m], forward_substitution(L_Phi[0:m, 0:m], B.T)))\
                              + backward_substitution(L_Phi.T[m:m+n+m, m:m+n+m],
                                                      forward_substitution(L_Phi[m:m+n+m, m:m+n+m], np.eye(n+m, n)))[0:n]
            # Y[i, i], i > 0
            elif i == j:
                Y[i*n:(i+1)*n, i*n:(i+1)*n] = (np.dot(np.hstack([A, B]), backward_substitution(L_Phi.T[m+(i-1)*(m+n):m+i*(m+n), m+(i-1)*(m+n):m+i*(m+n)],
                                                                    forward_substitution(L_Phi[m+(i-1)*(m+n):m+i*(m+n), m+(i-1)*(m+n):m+i*(m+n)], np.vstack([A.T, B.T]))))
                                               + backward_substitution(L_Phi.T[m+i*(m+n):m+i*(m+n)+n+m, m+i*(m+n):m+i*(m+n)+n+m],
                                                                        forward_substitution(L_Phi[m+i*(m+n):m+i*(m+n)+n+m, m+i*(m+n):m+i*(m+n)+n+m], np.eye(n+m, n)))[0:n])
            # Y[i, j] = Y[j, i]
            elif i != j and j <= T-1:
                Y[i*n:(i+1)*n, j*n:(j+1)*n] = (backward_substitution(L_Phi.T[m+i*(m+n):m+(i+1)*(m+n), m+i*(m+n):m+(i+1)*(m+n)],
                                                                    forward_substitution(L_Phi[m+i*(m+n):m+(i+1)*(m+n), m+i*(m+n):m+(i+1)*(m+n)], np.vstack([-A.T, -B.T]))))[0:5]
                Y[j*n:(j+1)*n].T[i*n:(i+1)*n] = Y[i*n:(i+1)*n].T[j*n:(j+1)*n].T
    return Y

def solve_lin_gs_structured(Phi, rd, rp, A, B, C, T, n, m):

    L_Phi = cholesky(Phi)

    Y = form_Y(L_Phi, A, B, T, n, m)
    beta = -rp + np.dot(C, backward_substitution(L_Phi.T, forward_substitution(L_Phi, rd)))

    # L_Y = np.zeros_like(Y)
    # L_Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = cholesky(Y[0*n:(0+1)*n].T[0*n:(0+1)*n])
    # for i in range(1, T):
    #
    #     L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n] = np.linalg.solve(L_Y[(i-1)*n:(i)*n].T[(i-1)*n:(i)*n], Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T)
    #     L_Y[i*n:(i+1)*n].T[i*n:(i+1)*n] = cholesky(Y[i*n:(i+1)*n].T[i*n:(i+1)*n] - np.dot(L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T, L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n]))

    # TODO einezelne Blöcke in L_Y berechnen sollte schneller gehen
    L_Y = cholesky(Y)

    delta_v = backward_substitution(L_Y.T, forward_substitution(L_Y, -beta))
    delta_z = backward_substitution(L_Phi.T, forward_substitution(L_Phi, -rd - np.dot(C.T, delta_v)))

    return np.vstack([delta_z, delta_v])

def forward_substitution(A, b):
    #TODO Prüfen ob A untere Dreiecksmatrix
    dim = [np.shape(A)[1], np.shape(b)[1]]
    x = np.array(np.eye(dim[0], dim[1]))
    for k in range(0, dim[0]):
        x[k] = (b[k] - np.dot(A[k, 0:k], x[0:k])) / A[k, k]
    return x

def backward_substitution(A, b):
    #TODO Prüfen ob A obere Dreiecksmatrix
    dim = [np.shape(A)[1], np.shape(b)[1]]
    x = np.array(np.eye(dim[0], dim[1]))
    for k in range(0, dim[0]):
        kh = dim[0] - k - 1
        x[kh] = (b[kh] - np.dot(A[kh, kh+1:], x[kh+1:])) / A[kh, kh]
    return x

def gradient(function, point, args=(), schritt=0.001):
    dim = np.shape(point)[0]
    grad = np.zeros([dim, 1])
    for i in range(dim):
        grad[i] = (function(point+(np.eye(dim)*schritt)[0:dim, [i]], *args)-function(point, *args))/schritt
    return grad

def backtracking_line_search(function, point, dir, args=()):
    # backtracking line search nach: TODO wonach?
    f_x = function(point, *args)
    grad_f = gradient(function, point, args, 0.000001)
    alpha = 0.4
    beta = 0.6
    st = 1
    while (np.isnan(function(point + st*dir, *args)) or
        function(point + st*dir, *args) > f_x + alpha*st*np.dot(grad_f.T, dir)):
        st = beta*st
        # print(st)
    return st