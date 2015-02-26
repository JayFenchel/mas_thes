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
    x=np.array(np.eye(npe, 1))
    x = backward_substitution(Rtilde, ctilde)
    return x

def solve_lin_gs_structured(Phi_not_inv, r, A, B, C, T, m, n, v):


    L_Phi = np.linalg.cholesky(Phi_not_inv)
############
    Phi = np.linalg.inv(Phi_not_inv)
############
    # TODO flogende Zeilen rausnehmen
    Y = np.eye(T*n, T*n)*0
    Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = np.dot(B, np.linalg.solve(L_Phi.T[0:m].T[0:m], np.linalg.solve(L_Phi.T[0:m].T[0:m].T, B.T))).T\
                                                + Phi[m:m+n].T[m:m+n].T
    for i in range(0, T):
        for j in range(i, i+2):
            #Y[0, 0]
            if i == 0 and j == 0:
                Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = np.dot(np.dot(B, Phi[0][0]), B.T).T\
                                                + Phi[m:m+n].T[m:m+n].T
            #Y[i, i]
            elif i == j:
                Y[i*n:(i+1)*n].T[i*n:(i+1)*n] = (np.dot(np.dot(A, Phi[m+(i-1)*(m+n):m+(i-1)*(m+n)+n].T[m+(i-1)*(m+n):m+(i-1)*(m+n)+n].T), A.T)
                                                + np.dot(np.dot(B, Phi[m+(i-1)*(m+n)+n:m+(i-1)*(m+n)+n+m].T[m+(i-1)*(m+n)+n:m+(i-1)*(m+n)+n+m].T), B.T)
                                                + Phi[m+i*(m+n):m+i*(m+n)+n].T[m+i*(m+n):m+i*(m+n)+n].T).T
            elif i != j and j <= T-1:
                Y[i*n:(i+1)*n].T[j*n:(j+1)*n] = -(np.dot(Phi[m+i*(m+n):m+i*(m+n)+n].T[m+i*(m+n):m+i*(m+n)+n].T, A.T)).T
                Y[j*n:(j+1)*n].T[i*n:(i+1)*n] = Y[i*n:(i+1)*n].T[j*n:(j+1)*n]

    beta = -r[T*(n+m):] + np.dot(np.dot(C, Phi), r[0:T*(n+m)])

    L_Y = np.zeros_like(Y)
    L_Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = cholesky(Y[0*n:(0+1)*n].T[0*n:(0+1)*n])
    for i in range(1, T):

        L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n] = np.linalg.solve(L_Y[(i-1)*n:(i)*n].T[(i-1)*n:(i)*n], Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T)
        L_Y[i*n:(i+1)*n].T[i*n:(i+1)*n] = cholesky(Y[i*n:(i+1)*n].T[i*n:(i+1)*n] - np.dot(L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T, L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n]))


    # L_Y = cholesky(Y)
    hilf=np.array(np.eye(np.shape(beta)[0],1))
    for k in range(0, 1, np.shape(beta)[0]):
        hilf[k] = (-beta[k] - np.dot(L_Y[k, 0:k], hilf[0 : k]))/L_Y[k, k]
    delta_v = np.linalg.solve(Y, -beta)
    for k in range(np.shape(beta)[0]-1, 0-1, -1):
        delta_v[k] = (hilf[k] - np.dot(L_Y.T[k, k+1:], hilf[k+1:]))/L_Y.T[k, k]
    # print delta_v[:10]
    delta_v = np.linalg.solve(Y, -beta)
    # print delta_v[:10]

    delta_z = np.linalg.solve(Phi_not_inv, -r[0:T*(n+m)] - np.dot(C.T, v))
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