__author__ = 'jayf'

from MyMath import forward_substitution
from MyMath import backward_substitution
from QuadraticProgram import QuadraticProgram
import numpy as np
from numpy.testing import assert_allclose
import unittest



class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.A = np.array([[  0.23996015,   0., 0.,   0., 0.],
                      [ -0.37221757,   1., 0.,   0., 0.],
                      [ -0.99008755,   0., 0.13885973,   0., 0.],
                      [-48.93540655, 64.1, 2.39923411,   1., 0.],
                      [  0.,           0., 0.,           0., 1.]])

        self.b = np.array([[-1.2346445 ,1],
                           [-1.43828223,1],
                           [-4.48282454,1],
                           [-1.79989043,1],
                           [1.,1]])

        T, n, m = 3, 5, 2
        self.test_qp = QuadraticProgram(T, n, m)

    # Test functions of MyMath
    def test_forward_substitution(self):

        x_f = forward_substitution(self.A, self.b)
        self.assertTrue((abs(np.dot(self.A, x_f)-self.b)).sum() < 1e-10,
                        'forward_substitution failed')

    def test_backward_substitution(self):

        x_b = backward_substitution(self.A.T, self.b)
        self.assertTrue((abs(np.dot(self.A.T, x_b)-self.b)).sum() < 1e-10,
                        'backward_substitution failed')

    # Test functions of Quadratic Program
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


if __name__ == '__main__':
    unittest.main()