#!/usr/bin/env python3
__author__ = 'jayf'

from ModelPredictiveControl.MyMath import solve_lin_gs
from ModelPredictiveControl.MyMath import form_Y
from ModelPredictiveControl.MyMath import solve_lin_gs_structured
from ModelPredictiveControl.MyMath import forward_substitution
from ModelPredictiveControl.MyMath import backward_substitution
from ModelPredictiveControl.MyMath import cholesky
from ModelPredictiveControl.SOCP import SOCP
import numpy as np
from numpy.testing import assert_allclose
import unittest


class MyTestCase(unittest.TestCase):

    def setUp(self):

        # For testing cholesky
        self.sym_matrix1 = np.array([[80., 2., 3., 4., 1., 1., 0.],
                                     [2., 90., 1., 6., 5., 4., 3.],
                                     [3., 1., 70., 1., 1., 1., 2.],
                                     [4., 6., 1., 40., 1., 1., 1.],
                                     [1., 5., 1., 1., 50., 5., 1.],
                                     [1., 4., 1., 1., 5., 60., 1.],
                                     [0., 3., 2., 1., 1., 1., 40.]])

        # For testing solve methods
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

        # For Tests requiring constraints
        self.Fx = np.array([[1, 2, 3, 4, 5],
                            [6, 7, 8, 9, 0],
                            [0, 0, 0, 1, 0]])
        self.fx = np.array([[1],
                            [2],
                            [3]])

        self.Fu = np.array([[2, 1],
                            [4, 3],
                            [6, 5]])
        self.fu = np.array([[0],
                            [1],
                            [2]])

        self.Ff = 0.5*self.Fx
        self.ff = 0.5*self.fx

        T, n, m = 3, 5, 2
        self.test_qp = SOCP(T, n, m)

    # Test functions of MyMath
    # TODO geeignete asserts für Toleranzbereiche wählen
    def test_solve_lin_gs_structured(self):

        self.test_qp.set_sys_dynamics(self.A, self.b)
        self.test_qp.set_lin_constraints(self.Fu, self.fu, self.Fx, self.fx,
                                     self.Ff, self.ff)
        x_test = np.array([[5, 1, 4, 2, 3]]).T
        z_test = np.array([[9, 0, 8, 1, 7, 2, 6, 3, 5, 4, 0, 9, 0, 8, 1, 7, 2,
                            6, 3, 5, 4]]).T
        v_test = np.array([[0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1]]).T
        zv_test = np.vstack([z_test, v_test])
        d = self.test_qp.form_d(x_test, zv_test)
        Phi = self.test_qp.form_Phi(d, z_test)
        rd, rp = self.test_qp.residual(x_test, zv_test)
        A, B, C = self.test_qp.A, self.test_qp.B, self.test_qp.C
        T, n, m = self.test_qp.T, self.test_qp.n, self.test_qp.m
        x = solve_lin_gs_structured(Phi, rd, rp, A, B, C, T, n, m)
        SS = np.hstack([np.vstack([Phi, C]), np.vstack([C.T, np.zeros([C.shape[0], C.shape[0]])])])
        self.assertTrue((abs(np.dot(SS, x)+np.vstack([rd, rp]))).sum() < 1e-10,
                        'solve_lin_gs_structured failed')

    # def test_form_Y(self):
    #     test_Phi_inv = np.array([[1., 0., 0., 0., 0., 0.],
    #                              [0., 2., 1., 0., 0., 0.],
    #                              [0., 4., 3., 0., 0., 0.],
    #                              [0., 0., 0., 2., 0., 0.],
    #                              [0., 0., 0., 0., 1., .5],
    #                              [0., 0., 0., 0., 3., 2.]])
    #
    #
    #     test_A = np.array([[0.5, 1.],
    #                        [2., 1.5]])
    #     test_B = np.array([[7.],
    #                        [6.]])
    #
    #     ref_Y = np.array([[51, 43],
    #                       [46, 39]])
    #
    #     Y = form_Y(test_Phi_inv, test_A, test_B, 2, 2, 1)
    #     print(Y)
    #     pass

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

    # TODO def test_gradient(self):
    #     def Quadrat(x):
    #         return (x*x).sum()

    # Test functions of QCQP
    # TODO def test_set_sys_dynamics(self):
    # TODO def test_set_weighting(self):
    # TODO def test_residual_norm(self):

    def test_set_ref_trajectory(self):
        x_ref_1 = np.array([[1],
                            [2],
                            [3],
                            [4],
                            [5]])
        x_ref_2 = np.array([[1, 2, 3, 4, 5],
                            [2, 2, 3, 4, 5],
                            [3, 3, 4, 4, 5],
                            [4, 4, 4, 4, 5],
                            [5, 5, 5, 5, 5]])

        ref_z_ref_1 = np.array([[0], [0], [1], [2], [3], [4], [5],
                                [0], [0], [1], [2], [3], [4], [5],
                                [0], [0], [1], [2], [3], [4], [5]])
        ref_update_1 = np.array([[1, 1]])
        ref_z_ref_2 = np.array([[0], [0], [1], [2], [3], [4], [5],
                                [0], [0], [2], [2], [3], [4], [5],
                                [0], [0], [3], [3], [4], [4], [5]])
        ref_update_2 = np.array([[4, 4, 4, 4, 5],
                                 [5, 5, 5, 5, 5]])
        self.test_qp.set_ref_trajectory(x_ref_1)
        self.assertTrue((self.test_qp.z_ref == ref_z_ref_1).all())
        self.assertTrue((self.test_qp.ref_update == ref_update_1).all())


        self.test_qp.set_ref_trajectory(x_ref_2)
        self.assertTrue((self.test_qp.z_ref == ref_z_ref_2).all())
        self.assertTrue((self.test_qp.ref_update == ref_update_2).all())

        # print(self.test_qp.z_ref)
        self.test_qp.update_ref_trajectory()
        # print(self.test_qp.z_ref)
        self.test_qp.update_ref_trajectory()
        # print(self.test_qp.z_ref)
        self.test_qp.update_ref_trajectory()
        # print(self.test_qp.z_ref)


if __name__ == '__main__':
    unittest.main()
