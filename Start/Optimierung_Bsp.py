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
#Startwerte für x(2), y(3), lambda(3)
x0 = np.matrix([0.5, 0.5]).T
y0 = [1., 1., 0.]
Y0 = np.diag(y0)
y0 = np.matrix(y0).T
lambda0 = [1., 1., 1.]
Lambda0 = np.diag(lambda0)
lambda0 = np.matrix(lambda0).T


x = x0
y = y0
lambda_ = lambda0
Y = Y0
Lambda = Lambda0
#alpha  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
for i in range(0, 10):
    my=y.T*lambda_/3.
    
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