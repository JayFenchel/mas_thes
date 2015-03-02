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

def form_Y(Phi_inv, A, B, T, n, m):

    # TODO find ungenauigkeit, irgendwo ist vllt noch ein kleiner fehler
    Y = np.zeros([T*n, T*n])

    for i in range(0, T):
        for j in range(i, i+2):
            # Y[0, 0]
            if i == 0 and j == 0:
                Y[0:n, 0:n] = np.dot(np.dot(B, Phi_inv[0:m, 0:m]), B.T)\
                                                + Phi_inv[m:m+n, m:m+n]
            # Y[i, i], i > 0
            elif i == j:
                # print(Phi_inv[m+(i-1)*(m+n):m+(i-1)*(m+n)+n, m+(i-1)*(m+n)+n:m+i*(m+n)])
                Y[i*n:(i+1)*n, i*n:(i+1)*n] = (np.dot(np.dot(A, Phi_inv[m+(i-1)*(m+n):m+(i-1)*(m+n)+n, m+(i-1)*(m+n):m+(i-1)*(m+n)+n]), A.T)
                                                + np.dot(A, np.dot(Phi_inv[m+(i-1)*(m+n):m+(i-1)*(m+n)+n, m+(i-1)*(m+n)+n:m+i*(m+n)], B.T))
                                                + np.dot(B, np.dot(Phi_inv[m+(i-1)*(m+n):m+(i-1)*(m+n)+n, m+(i-1)*(m+n)+n:m+i*(m+n)].T, A.T))
                                                + np.dot(np.dot(B, Phi_inv[m+(i-1)*(m+n)+n:m+i*(m+n), m+(i-1)*(m+n)+n:m+i*(m+n)]), B.T)
                                                + Phi_inv[m+i*(m+n):m+i*(m+n)+n, m+i*(m+n):m+i*(m+n)+n])
            elif i != j and j <= T-1:
                Y[i*n:(i+1)*n, j*n:(j+1)*n] = (-(np.dot(Phi_inv[m+i*(m+n):m+i*(m+n)+n, m+i*(m+n):m+i*(m+n)+n], A.T))
                                               -(np.dot(Phi_inv[m+i*(m+n):m+i*(m+n)+n, m+i*(m+n)+n:m+(i+1)*(m+n)], B.T)))
                Y[j*n:(j+1)*n].T[i*n:(i+1)*n] = Y[i*n:(i+1)*n].T[j*n:(j+1)*n].T
    # print(Y)
    return Y

def solve_lin_gs_with_Y(Phi, C, rd, rp):

    Phi_inv = np.linalg.inv(Phi)

    Y = np.dot(C, np.dot(Phi_inv, C.T))
    beta = - rp + np.dot(C, np.dot(Phi_inv, rd))

    delta_v = np.linalg.solve(Y, -beta)
    delta_z = np.linalg.solve(Phi, -rd-np.dot(C.T, delta_v))

    x = np.vstack([delta_z, delta_v])
    return x

def solve_lin_gs_structured(Phi, rd, rp, A, B, C, T, n, m):

    Phi_inv = np.linalg.inv(Phi)
    L_Phi = cholesky(Phi)
    Y = form_Y(Phi_inv, A, B, T, n, m)

    beta = -rp + np.dot(np.dot(C, Phi_inv), rd)

    # L_Y = np.zeros_like(Y)
    # L_Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = cholesky(Y[0*n:(0+1)*n].T[0*n:(0+1)*n])
    # for i in range(1, T):
    #
    #     L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n] = np.linalg.solve(L_Y[(i-1)*n:(i)*n].T[(i-1)*n:(i)*n], Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T)
    #     L_Y[i*n:(i+1)*n].T[i*n:(i+1)*n] = cholesky(Y[i*n:(i+1)*n].T[i*n:(i+1)*n] - np.dot(L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T, L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n]))

    # TODO einezelne Blöcke in L_Y berechnen sollte schneller gehen
    L_Y = cholesky(Y)
    delta_v_hilf = forward_substitution(L_Y, -beta)
    delta_v = backward_substitution(L_Y.T, delta_v_hilf)

    delta_z_hilf = forward_substitution(L_Phi, -rd - np.dot(C.T, delta_v))
    delta_z = backward_substitution(L_Phi.T, delta_z_hilf)

    x = np.vstack([delta_z, delta_v])
    return x

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