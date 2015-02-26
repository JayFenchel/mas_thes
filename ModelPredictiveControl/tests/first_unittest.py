#!/usr/bin/env python3
__author__ = 'jayf'

from ModelPredictiveControl.MyMath import solve_lin_gs
from ModelPredictiveControl.MyMath import solve_lin_gs_structured
from ModelPredictiveControl.MyMath import forward_substitution
from ModelPredictiveControl.MyMath import backward_substitution
from ModelPredictiveControl.MyMath import cholesky
from ModelPredictiveControl.QuadraticProgram import QuadraticProgram
import numpy as np
from numpy.testing import assert_allclose
import unittest


class MyTestCase(unittest.TestCase):

    def setUp(self):

        self.sym_matrix1 = np.array([[80., 2., 3., 4., 1., 1., 0.],
                                     [2., 90., 1., 6., 5., 4., 3.],
                                     [3., 1., 70., 1., 1., 1., 2.],
                                     [4., 6., 1., 40., 1., 1., 1.],
                                     [1., 5., 1., 1., 50., 5., 1.],
                                     [1., 4., 1., 1., 5., 60., 1.],
                                     [0., 3., 2., 1., 1., 1., 40.]])

        self.A = np.array([[0.23996015,    0.,  0.,           0., 0.],
                           [-0.37221757,   1.,  0.,           0., 0.],
                           [-0.99008755,   0.,  0.13885973,   0., 0.],
                           [-48.93540655, 64.1, 2.39923411,   1., 0.],
                           [0.,            0.,  0.,           0., 1.]])

        self.b = np.array([[-1.2346445,  1.],
                           [-1.43828223, 1.],
                           [-4.48282454, 1.],
                           [-1.79989043, 1.],
                           [1.,          1.]])

        T, n, m = 3, 5, 2
        self.test_qp = QuadraticProgram(T, n, m)

    # Test functions of MyMath
    # TODO geeignete asserts für Toleranzbereiche wählen
    # def test_solve_lin_gs_structured(self):
    #     x = solve_lin_gs_structured(self.A, self.b)
    #     self.assertTrue((abs(np.dot(self.A, x)-self.b)).sum() < 1e-10,
    #                     'solve_lin_gs_structured failed')

    def test_solve_lin_gs(self):
        x = solve_lin_gs(self.A, self.b)
        self.assertTrue((abs(np.dot(self.A, x)-self.b)).sum() < 1e-10,
                        'solve_lin_gs failed')

    def test_cholesky(self):
        ref = np.linalg.cholesky(self.sym_matrix1)
        c = cholesky(self.sym_matrix1)
        self.assertTrue((abs(c - ref)).sum() < 1e-10,
                        'cholesky decomposition failed')
        c = cholesky(self.sym_matrix1)
        self.assertTrue((abs(c - ref)).sum() < 1e-10,
                        'cholesky decomposition failed (second)')

    def test_forward_substitution(self):
        x_f = forward_substitution(self.A, self.b)
        self.assertTrue((abs(np.dot(self.A, x_f)-self.b)).sum() < 1e-10,
                        'forward_substitution failed')

    def test_backward_substitution(self):
        x_b = backward_substitution(self.A.T, self.b)
        self.assertTrue((abs(np.dot(self.A.T, x_b)-self.b)).sum() < 1e-10,
                        'backward_substitution failed')

    # Test functions of Quadratic Program
    # TODO def test_set_sys_dynamics(self):
    # TODO def test_set_weighting(self):
    def test_set_constraints(self):

        Fx = np.array([[1, 2, 3, 4, 5],
                       [6, 7, 8, 9, 0],
                       [0, 0, 0, 1, 0]])
        fx = np.array([[1],
                       [2],
                       [3]])

        Fu = np.array([[2, 1],
                       [4, 3],
                       [6, 5]])
        fu = np.array([[0],
                       [1],
                       [2]])

        Ff = 0.5*Fx
        ff = 0.5*fx

        self.test_qp.set_constraints(Fu, fu, Fx, fx, Ff, ff)
        ref = np.array(
            [[2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 1, 2, 3, 4, 5, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 6, 7, 8, 9, 0, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 1, 0, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 2, 1, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 8, 9, 0, 4, 3, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 6, 5, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1.0, 1.5, 2.0, 2.5],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.0, 3.5, 4.0, 4.5,   0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,   0,   0, 0.5,   0]])

        self.assertTrue((self.test_qp.P == ref).all(), 'False P-matrix')
        # TODO Test h-vector

        # TODO Test form_d
        z_test = np.ones([21, 1])
        print(z_test)

if __name__ == '__main__':
    unittest.main()
