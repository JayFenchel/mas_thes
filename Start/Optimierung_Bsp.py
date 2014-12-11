#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:46:04 2014

@author: jayf
"""

from cvxopt import matrix
import numpy as np

def function_parameter():
    
    Q = 2*matrix([ [2, .5], [.5, 1] ])  #System
    p = matrix([1.0, 1.0])              #System
    G = matrix([[-1.0,0.0],[0.0,-1.0]]) #UN
    h = matrix([0.0,0.0])               #UN
    A = matrix([1.0, 1.0], (1,2))       #GN
    b = matrix(1.0)                     #GN

    # Notation Buch
    G_ip = np.matrix(Q)
    c_ip = np.matrix(p)
    A_ip = np.vstack([-np.matrix(G), np.matrix(A)])
    b_ip = np.vstack([np.matrix(h), np.matrix(b)])
    
    return G_ip, c_ip, A_ip, b_ip


G, c, A, b = function_parameter()
m = np.shape(A)[0]
#Startwerte für x(2), y(3), lambda(3)
#Compute (x0, y0, lambda_0) with (y0,lambda_0) >0
x0 = np.matrix([0.5, 0.5]).T
y0 = np.matrix([1., 1., 0.]).T
lambda_0 = np.matrix([1., 1., 1.]).T

#zum besseren Anwenden des Algorithmus
xk = x0
yk = y0
lambda_k = lambda_0

for k in range(0, 10):
    
    #Set (x, y, lambda_) = ((xk, yk, lambda_k))
    x, y, lambda_ = xk, yk, lambda_k
    
    #Calculate my = y.T*lambda_/m
    my = y.T*lambda_/m
    
    #Solve (...) with sigma = 0 for (delta_x_aff, delta_y_aff, delta_lambda__aff)
    #
    #
    #
    #...........................
    
    V = np.hstack([np.vstack([G, A, np.zeros_like(A)]),
                   np.vstack([np.zeros_like(A.T), -np.matrix(np.eye(3)), Lambda]),
                   np.vstack([-A.T, np.zeros_like(Y), Y])])
                   
    H = np.vstack([-(G*x - A.T*lambda_ + c),
                   -(A*x - y - b),
                   -Lambda*Y*np.matrix([1., 1., 0. ]).T + 0*my.flat[0]*np.matrix([1., 1., 1. ]).T])
                   
    #lin GS lösen
    delta = np.linalg.solve(V, H)
    
    #und (x+, y+, lambda+) = (x, y, lambda) + delta(x, y, lambda)
    
    alpha = 0.9
    
    x += alpha*delta[0:2]
    y += alpha*delta[2:5]
    lambda_ += alpha*delta[5:8]
    
    print(np.vstack([y[0:2],lambda_]).min())
#    while np.vstack([y[0:2],lambda_]).min() <= 0:
#        print(np.vstack([y[0:2],lambda_]).min())
#        alpha *= 0.9
#        x -= alpha*delta[0:2]
#        y -= alpha*delta[2:5]
#        lambda_ -= alpha*delta[5:8]
#        
    Y = np.diag([y.flat[0],y.flat[1],y.flat[2]])
    Lambda = np.diag([lambda_.flat[0],lambda_.flat[1],lambda_.flat[2]])
    print(alpha)
    
    iter_ = np.vstack([x, y, lambda_])

print(delta)
print(iter_)
input()

def matrix_diag(a):
    
    return np.diag([a.flat[0],a.flat[1],a.flat[2]])