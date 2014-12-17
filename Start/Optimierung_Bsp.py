#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:46:04 2014

@author: jayf
"""

#~ from cvxopt import matrix
import numpy as np
from pprint import pprint

class GS:
    
    def __init__(self, get_gleichung, get_G, get_c, get_A, get_b):
        G = get_G
        c = get_c
        A = get_A
        b = get_b
        if get_gleichung == 1 or get_gleichung == 2:
            self.gleichung = get_gleichung
        else:
            self.gleichung = 1
        
        
    def solve(self, x, y, lambda_, sigma, my, delta_aff):
        
        Y = matrix_diag(y)
        Lambda = matrix_diag(lambda_)

        V = np.hstack([
                np.vstack([G, A, np.zeros_like(A)]),
                np.vstack([np.zeros_like(A.T), -np.array(np.eye(3)), Lambda]),
                np.vstack([-A.T, np.zeros_like(Y), Y])])
        
        if self.gleichung == 1:           
            H = np.vstack([
                    -(np.dot(G, x) - np.dot(A.T, lambda_) + c),
                    -(np.dot(A, x) - y - b),
                    -np.dot(np.dot(Lambda, Y), np.array([[1.], [1.], [0.] ])) 
                        + sigma*my*np.array([[1.], [1.], [0.] ])])
                           
        if self.gleichung == 2:
            delta_aff_Y = matrix_diag(delta_aff[2:5])
            delta_aff_Lambda = matrix_diag(delta_aff[5:8])
            
            H = np.vstack([
                    -(np.dot(G, x) - np.dot(A.T, lambda_) + c),
                    -(np.dot(A, x) - y - b),
                    -np.dot(np.dot(Lambda, Y), np.array([[1.], [1.], [0.] ])) 
                        - np.dot(np.dot(delta_aff_Lambda, delta_aff_Y), np.array([[1.], [1.], [0.] ])) 
                        + sigma*my*np.array([[1.], [1.], [0.] ])])
        
        return np.linalg.solve(V, H)    # lin GS lösen

def matrix_diag(a):
    
    return np.diag(a.T[0])

def function_parameter():
    
    #~ Q = 2*matrix([ [2, .5], [.5, 1] ])   # KF
    #~ p = matrix([1.0, 1.0])               # KF
    #~ G = matrix([[-1.0, 0.0],[0.0,-1.0]]) # UN
    #~ h = matrix([0.0, 0.0])               # UN
    #~ A = matrix([1.0, 1.0], (1,2))        # GN
    #~ b = matrix(1.0)                      # GN
    
    Q = [ [4., 1.], [1., 2.] ]              # KF
    p = [[1.0], [1.0]]                      # KF
    G = [[-1.0,0.0],[0.0,-1.0]]             # UN
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
m = np.shape(A)[0] # Number of Zeilen
# Startwerte für x(2), y(3), lambda(3)
# Compute (x0, y0, lambda_0) with (y0,lambda_0) >0
x0 = np.array([[0.5], [0.5]])
y0 = np.array([[1.], [1.], [0.]])
lambda_0 = np.array([[1.], [1.], [1.]])

# zum besseren Anwenden des Algorithmus
xk = x0
yk = y0
lambda_k = lambda_0

GS1 = GS(1, G, c, A, b)
GS2 = GS(2, G, c, A, b)

for k in range(0, 10):
    
    # Set (x, y, lambda_) = ((xk, yk, lambda_k))
    x, y, lambda_ = xk, yk, lambda_k
    
    # Calculate my = y.T*lambda_/m
    my = np.dot(y.T, lambda_)/m
    
    # Solve (...) with sigma = 0 for (delta_x_aff, delta_y_aff,
    # delta_lambda__aff)
    sigma = 0
    delta_aff = GS1.solve(x, y, lambda_, sigma, my, 0)
    
    # Calculate alpha_aff_Dach = max(...)
    alpha_aff_dach = 0.9 # TODO
    
    # Calculate my_dach = ...
    my_aff = (np.dot((y + alpha_aff_dach*delta_aff[2:5]).T, 
             (lambda_ + alpha_aff_dach*delta_aff[5:8])) / m)

    # Set centering parameter to sigma = ...
    sigma = np.power(my_aff/my, 3)
    
    # Solve (...) for (delta_x, delta_y, delta_lambda_)
    delta = GS2.solve(x, y, lambda_, sigma, my, delta_aff)
    
    # Choose tau_k ... and set alpha_dach = ... (see(16.66))    
    alpha_dach = 0.9 # TODO
    
    # Set (x+, y+, lambda+) = (x, y, lambda) + alpha_dach
    # * delta(x, y, lambda)    
    xk = x + alpha_dach*delta[0:2]
    yk = y + alpha_dach*delta[2:5]
    lambda_k = lambda_ + alpha_dach*delta[5:8]
    
    print(np.vstack([x, y, lambda_]).T)

input()

