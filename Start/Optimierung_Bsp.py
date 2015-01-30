#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:46:04 2014

@author: jayf
"""

import numpy as np


class LinearSystem:

    def __init__(self, get_gleichung, get_G, get_c, get_A, get_b):
        G = get_G
        c = get_c
        A = get_A
        b = get_b
        if get_gleichung == 1 or get_gleichung == 2:
            self.gleichung = get_gleichung
        else:
            self.gleichung = 1

    def solve_lin_gs(self, A, b):

        neta, npe = A.shape  # neta - # of Zeilen, npe - # of Spalten in A
        Q, R = householder(A)
        # hhhh = cholesky(A)
        c = np.dot(Q, b)
        Rtilde = R[0:npe, 0:npe]
        ctilde = c[0:npe]
        x=np.array(np.eye(npe, 1))
        for k in range(npe-1, 0-1, -1):
            x[k] = (ctilde[k] - np.dot(Rtilde[k, k+1:npe], x[k+1:npe]))/Rtilde[k, k]

        return x

    def solve(self, x, y, lambda_, sigma, my, delta_aff):
        
        Y = matrix_diag(y)

        Y_inv = np.linalg.pinv(Y)
        Y_inv[2][2] = 1
        Lambda = matrix_diag(lambda_)

        V = np.hstack([
                np.vstack([G, np.zeros_like(A), -A]),
                np.vstack([np.zeros_like(A.T), np.dot(Y_inv, Lambda), +np.array(np.eye(3))]),
                np.vstack([-A.T, np.dot(Y, Y_inv), np.zeros_like(Y)])])
        
        if self.gleichung == 1:
            H = np.vstack([
                    -(np.dot(G, x) - np.dot(A.T, lambda_) + c),
                    -np.dot(np.dot(np.dot(Y_inv, Lambda), Y), np.array([[1.], [1.], [0] ]))
                        + np.dot(Y_inv, sigma*my*np.array([[1.], [1.], [0] ])),
                    +(np.dot(A, x) - y - b)])

        if self.gleichung == 2:
            delta_aff_Y = matrix_diag(delta_aff[2:5])
            delta_aff_Lambda = matrix_diag(delta_aff[5:8])
            
            H = np.vstack([
                    -(np.dot(G, x) - np.dot(A.T, lambda_) + c),
                    -np.dot(np.dot(np.dot(Y_inv,Lambda), Y), np.array([[1.], [1.], [0] ]))
                        - np.dot(np.dot(np.dot(Y_inv, delta_aff_Lambda), delta_aff_Y), np.array([[1.], [1.], [0] ]))
                        + np.dot(Y_inv,sigma*my*np.array([[1.], [1.], [0] ])),
                    +(np.dot(A, x) - y - b)])

        # return np.linalg.solve(V, H)    # lin GS lösen
        return self.solve_lin_gs(V, H)    # lin GS lösen


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

def cholesky(a):
    # TODO
    n=a.shape[0]
    print(a)
    for i in range(0, n):
        for j in range(0, i):
            sum = a[i, j]
            for k in range(0, j-1):
                sum += -a[i, k]*a[j, k]
            if i>j:
                a[i, j] = sum/np.abs(a[j, j])
            else:
                if sum > 0:
                    a[i, i] = np.sqrt(sum)
                else:
                    print('ERROR')
    #
    #
    # for k in range(0, 8):
    #     for j in range(0, k-1):
    #         a[k, k] = a[k, k] - a[k, j]*a[k, j]
    #     a[k, k] = np.sqrt(np.abs(a[k, k]))
    #     for i in range(k+1, 8):
    #         for j in range(0, k-1):
    #             a[i, k] = a[i, k] - a[i, j]*a[k, j]
    #         a[i, k] = a[i, k]/a[k, k]

    # Vergleich mit Housholder Transformation
    G = None
    # print(a)
    return(G)


def vector_norm(a):
    n = 0
    for i in range(0, a.shape[0]):
        n += a[i]**2
    return np.sqrt(n)

def matrix_diag(a):
    
    return np.diag(a.T[0])


def function_parameter():
    
    # Q = 2*matrix([ [2, .5], [.5, 1] ])   # KF
    # p = matrix([1.0, 1.0])               # KF
    # G = matrix([[-1.0, 0.0],[0.0,-1.0]]) # UN
    # h = matrix([0.0, 0.0])               # UN
    # A = matrix([1.0, 1.0], (1,2))        # GN
    # b = matrix(1.0)                      # GN
    
    Q = [[4., 1.], [1., 2.]]              # KF
    p = [[1.0], [1.0]]                      # KF
    G = [[-1.0, 0.0], [0.0, -1.0]]             # UN
    h = [[0.0], [0.0]]                      # UN
    A = [[1.0, 1.0]]                        # GN
    b = 1.0                                 # GN

    # Notation Buch
    G_ip = np.array(Q)
    c_ip = np.array(p)
    A_ip = np.vstack([-np.array(G), np.array(A)])
    b_ip = np.vstack([np.array(h), np.array(b)])
    
    return G_ip, c_ip, A_ip, b_ip


G, c, A, b = function_parameter()
m = np.shape(A)[0]  # Number of Zeilen
# Startwerte für x(2), y(3), lambda(3)
# Compute (x0, y0, lambda_0) with (y0,lambda_0) >0
x0 = np.array([[0.5], [0.5]])
y0 = np.array([[1.], [1.], [0.]])
lambda_0 = np.array([[1.], [1.], [1.]])

# zum besseren Anwenden des Algorithmus
xk = x0
yk = y0
lambda_k = lambda_0

GS1 = LinearSystem(1, G, c, A, b)
GS2 = LinearSystem(2, G, c, A, b)

for k in range(0, 6):
    
    # Set (x, y, lambda_) = ((xk, yk, lambda_k))
    x, y, lambda_ = xk, yk, lambda_k
    
    # Calculate my = y.T*lambda_/m
    my = np.dot(y.T, lambda_)/m
    
    # Solve (...) with sigma = 0 for (delta_x_aff, delta_y_aff,
    # delta_lambda__aff)
    sigma = 0
    delta_aff = GS1.solve(x, y, lambda_, sigma, my, 0)
    
    # Calculate alpha_aff_Dach = max(...)
    alpha_test = 1
    for i in range(0, 10):
        y_test = y + alpha_test*delta_aff[2:5]
        lambda_test = lambda_ + alpha_test*delta_aff[5:8]
        if np.vstack([y_test[0:2], lambda_test]).min() <= 0:
            alpha_test += -0.1
        else:
            break

    alpha_aff_dach = alpha_test  # TODO
    # Calculate my_dach = ...
    my_aff = (np.dot((y + alpha_aff_dach*delta_aff[2:5]).T, 
             (lambda_ + alpha_aff_dach*delta_aff[5:8])) / m)

    # Set centering parameter to sigma = ...
    sigma = np.power(my_aff/my, 3)
    
    # Solve (...) for (delta_x, delta_y, delta_lambda_)
    delta = GS2.solve(x, y, lambda_, sigma, my, delta_aff)
    
    # Choose tau_k ... and set alpha_dach = ... (see(16.66))    
    tau_k = 1  # tau_k is element(0, 1)

    alpha_pri_test = 1
    for i in range(0, 10):
        y_test = y + alpha_pri_test*delta[2:5]
        if y_test[0:2].min() <= 0:
            alpha_pri_test += -0.1
        else:
            break

    alpha_dual_test = 1
    for i in range(0,10):
        lambda_test = lambda_ + alpha_dual_test*delta[5:8]
        # print(np.vstack([y[0:2],lambda_]))
        if lambda_test.min() <= 0:
            alpha_dual_test += -0.1
        else:
            break

    alpha_dach = min(alpha_pri_test, alpha_dual_test) # TODO
    # Set (x+, y+, lambda+) = (x, y, lambda) + alpha_dach
    # * delta(x, y, lambda)    
    xk = x + alpha_dach*delta[0:2]
    yk = y + alpha_dach*delta[2:5]
    lambda_k = lambda_ + alpha_dach*delta[5:8]

    epsilon = 1e-10
    if abs(yk[2]) >= epsilon:
        print('Gleichungsnebenbed. nicht erfüllt!')
    if np.vstack([yk[0:2], lambda_k]).min() <= 0:
        print('Größer-gleich-Null-Bed. nicht erfüllt')
    print(xk.T)


