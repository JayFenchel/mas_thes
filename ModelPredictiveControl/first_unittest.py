__author__ = 'jayf'

from QuadraticProgram import QuadraticProgram
import numpy as np
from numpy.testing import assert_allclose
import unittest



class MyTestCase(unittest.TestCase):

    def setUp(self):
        T, n, m = 3, 5, 2
        self.test_qp = QuadraticProgram(T, n, m)

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

        self.test_qp.set_constraints(Fu, fu, Fx, fx, Ff, ff)


        self.assertTrue((self.test_qp.P == ref).all())


if __name__ == '__main__':
    unittest.main()
