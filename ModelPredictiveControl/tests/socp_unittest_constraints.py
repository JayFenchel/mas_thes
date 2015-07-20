#!/usr/bin/env python3
__author__ = 'jayf'

from ModelPredictiveControl.SOCP import SOCP
import numpy as np
import unittest

class MyTestCase(unittest.TestCase):

    def setUp(self):
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

        self.ref_P = np.array(
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

        T, n, m = 3, 5, 2
        self.test_qp = SOCP(T, n, m)

    def test_set_lin_constraints(self):

        self.test_qp.set_lin_constraints(self.Fu, self.fu, self.Fx, self.fx,
                                         self.Ff, self.ff)

        self.assertTrue((self.test_qp.P_of_zk(None) == self.ref_P).all(), 'False P-matrix')

        ref_h = np.array([[1, 3, 5, 1, 3, 5, 1, 3, 5, 0.5, 1, 1.5]]).T
        self. assertTrue((self.test_qp.h_of_xk(np.array([[0], [0], [0], [0], [0]])) == ref_h).all(), 'False h-vector')

        x_test = np.array([[5, 1, 4, 2, 3]]).T
        z_test = np.array([[9, 0, 8, 1, 7, 2, 6, 3, 5, 4, 0, 9, 0, 8, 1, 7, 2,
                            6, 3, 5, 4]]).T
        v_test = np.array([[0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1]]).T
        zv_test = np.vstack([z_test, v_test])

        ref_d = np.array([[-0.01694915254, -0.00833333333, -0.01960784314,
                           -0.01265822785, -0.00653594771, -0.02500000000,
                           -0.01265822785, -0.00847457627, -0.02777777778,
                           -0.03225806452, -0.01652892562, -1.00000000000]]).T
        self.assertTrue((abs(self.test_qp.form_d(x_test, zv_test) - ref_d)).sum() < 1e-10,
                        'False d-vector')
        # TODO Test residual(), dazu noch C-Matrix setzen
        # # Mal angenommen, das alte residual was richtig
        # self.assertTrue((np.vstack(self.test_qp.old_residual(x_test, zv_test)) ==
        #                 np.vstack(self.test_qp.residual(x_test, zv_test))).all(),
        #                 'Residual changed')

    # def test_eval_of_non_lin_constraints(self):
    #
    #     self.test_qp.set_lin_constraints(self.Fu, self.fu, self.Fx, self.fx,
    #                                      self.Ff, self.ff)
    #
    #
    #     self.test_qp.add_qc(F_qc=None, alpha=None)
    #
    #
    #     socc_A_test =
    #     self.test_qp.add_socc(socc_A=socc_A_test, socc_c=M21.T,
    #             socc_b=np.zeros_like(M21.T), socc_d=c)
    #
    #     z_test = np.array([[9, 0, 8, 1, 7, 2, 6, 3, 5, 4, 0, 9, 0, 8, 1, 7, 2,
    #                         6, 3, 5, 4]]).T
    #     ref_P = np.zeros([np.shape(self.ref_P)[0]+self.test_qp.T-1, np.shape(self.ref_P)[1]])
    #     # TODO auch extra Zeile nach erstem Fu?
    #     ref_P[0:6] = self.ref_P[0:6]
    #     ref_p[6] =
    #     ref_P[7:10] = self.ref_P[6:9]
    #     ref_p[10] =
    #     ref_P[11:14] = self.ref_P[9:12]
    #     self.assertTrue((self.test_qp.P_of_zk(z_test) == ref_P).all(),
    #                     'P-matrix does not changed in the right way')
    #     print(ref_P)


if __name__ == '__main__':
    unittest.main()