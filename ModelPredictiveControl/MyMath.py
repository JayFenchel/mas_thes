#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'
# Build with: python3.4 setup.py build_ext --inplace

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

def form_Y(Phi, A, B, T, n, m):

    # TODO find ungenauigkeit, irgendwo ist vllt noch ein kleiner fehler
    Y = np.zeros([T*n, T*n])
    L_Y = np.zeros_like(Y)
    L_Phi = np.zeros_like(Phi)
    L_Phi[0:m, 0:m] = cholesky(Phi[0:m, 0:m])

    for i in range(0, T):
        L_Phi[m+i*(m+n):m+(i+1)*(m+n), m+i*(m+n):m+(i+1)*(m+n)] =\
                cholesky(Phi[m+i*(m+n):m+(i+1)*(m+n), m+i*(m+n):m+(i+1)*(m+n)])
        bl_i = L_Phi[m+i*(m+n):m+(i+1)*(m+n), m+i*(m+n):m+(i+1)*(m+n)]
        for j in range(i, i+2):

            # Y[0, 0]
            if i == 0 and j == 0:
                r0 = L_Phi[0:m, 0:m]
                Y[0:n, 0:n] = np.dot(B, backward_substitution(r0.T, forward_substitution(r0, B.T)))\
                              + backward_substitution(bl_i.T, forward_substitution(bl_i, np.eye(n+m, n)))[0:n]
                L_Y[0:n, 0:n] = cholesky(Y[0:n, 0:n])
            # Y[i, i], i > 0
            elif i == j:
                bl_i_m1 = L_Phi[m+(i-1)*(m+n):m+i*(m+n), m+(i-1)*(m+n):m+i*(m+n)]
                Y[i*n:(i+1)*n, i*n:(i+1)*n] = (np.dot(np.hstack([A, B]), backward_substitution(bl_i_m1.T,
                                                                    forward_substitution(bl_i_m1, np.vstack([A.T, B.T]))))
                                               + backward_substitution(bl_i.T,
                                                                        forward_substitution(bl_i, np.eye(n+m, n)))[0:n])
            # Y[i, j] = Y[j, i]
            elif i != j and j <= T-1:
                Y[i*n:(i+1)*n, j*n:(j+1)*n] = (backward_substitution(bl_i.T,
                                                                    forward_substitution(bl_i, np.vstack([-A.T, -B.T]))))[0:n]
                Y[j*n:(j+1)*n].T[i*n:(i+1)*n] = Y[i*n:(i+1)*n].T[j*n:(j+1)*n].T
    return Y, L_Phi

def solve_lin_gs_structured(Phi, rd, rp, A, B, C, T, n, m, reg=0):
    # reg=epsilon_for_regularization, reg=0 => no regularization

    # L_Phi = cholesky(Phi)
    # regularization (+epsilon*I)
    Phi += reg*np.eye(np.shape(Phi)[0])  # Phi has to be quadratic
    Y, L_Phi = form_Y(Phi, A, B, T, n, m)
    beta = -rp + np.dot(C, backward_substitution(L_Phi.T, forward_substitution(L_Phi, rd)))

    # L_Y = np.zeros_like(Y)
    # L_Y[0*n:(0+1)*n].T[0*n:(0+1)*n] = cholesky(Y[0*n:(0+1)*n].T[0*n:(0+1)*n])
    # for i in range(1, T):
    #
    #     L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n] = np.linalg.solve(L_Y[(i-1)*n:(i)*n].T[(i-1)*n:(i)*n], Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T)
    #     L_Y[i*n:(i+1)*n].T[i*n:(i+1)*n] = cholesky(Y[i*n:(i+1)*n].T[i*n:(i+1)*n] - np.dot(L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n].T, L_Y[i*n:(i+1)*n].T[(i-1)*n:(i)*n]))

    # regularization (-epsilon*I)
    Y += reg*np.eye(np.shape(Y)[0])
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

def gradient(function, point, args=(), step=0.001):
    dim = np.shape(point)[0]
    steps = np.eye(dim)*step
    f_x = function(point, *args)
    grad = np.zeros([dim, 1])
    for i in range(dim):
        step_i = steps[0:dim, [i]]
        # if point+step_i is not feasible
        if np.isnan(function(point+step_i, *args)):
            grad[i] = (f_x - function(point-step_i, *args))/step
        else:
            grad[i] = (function(point+step_i, *args) - f_x)/step
    return grad

# TODO Anpassungen in Quadratic Programm, damit gradient_better wieder genutzt
# werden kann
def gradient_better(function, point, args=(), schritt=0.001):
    dim = np.shape(point)[0]
    grad = np.zeros([dim, 1])
    grad[:, 0] = (function(point+(np.eye(dim)*schritt)[0:dim, :], *args) - function(point, *args))/schritt
    return grad

# def hessian(args=(), schritt=0.001):
#     vec1 = np.array([[0.1, 0.], [0., 0.1]])
#     vec2 = np.array([[[0.1], [0.]], [[0.], [0.1]]])
#     print(vec1+vec2)
#     print('stop')
#     dim = np.shape(point)[0]
#     grad = np.zeros([dim, 1])
#     function_values = function(point+(np.eye(dim)*schritt)[0:dim, :], *args)
#
#     grad[:, 0] = (function_values - function(point, *args))/schritt
#     hess = 0
#     return hess

def backtracking_line_search(function, point, drect, args=(), step = 1e-08):
    # backtracking line search nach:
    # Boyd, Convex Optimization, Seite 464, Algorithm 9.2
    f_x = function(point, *args)
    grad_f = gradient(function, point, args, step)
    alpha = 0.4
    beta = 0.6
    st = 1
    grad_in_dir = np.dot(grad_f.T, drect)
    print('gradient in search direction = %s' % grad_in_dir)
    if np.isnan(grad_in_dir):
        print('Gradient konnte nicht berechnet werden, da delta in den unzulässigen Bereich geht.')
        exit()
    if grad_in_dir > 0:
        # Falls Kostenfunktion in Schrittrichtung am aktuellen Punkt ansteigend
        print('Gradient is positive, no improvement possible.')
        alpha = 0  # Funktionswert muss besser werden als aktuell
    f_x_step = function(point + st*drect, *args)
    while np.isnan(f_x_step) or f_x_step > f_x + alpha*st*grad_in_dir:
        st *= beta
        f_x_step = function(point + st*drect, *args)
        # Wenn Schrittweite immer kleiner, wahrscheinlich wegen nummerischem
        # Fehler Werte in Bedingung gleich groß für alpha = 0 Fall
        # print(st)
    return st

def backtracking_line_search_quick_and_dirty(function, point, drect, args=(), step = 0.000001):
    # only use gradient in search direction
    f_x = function(point, *args)
    grad_in_dir = function(point+step*drect, *args) - f_x
    alpha = 0.4
    beta = 0.6
    st = 1
    print(grad_in_dir)
    while (np.isnan(function(point + st*drect, *args)) or
            function(point + st*drect, *args) > f_x + alpha*st*grad_in_dir):
        st = beta*st
    return st

# hessian()