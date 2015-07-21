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

        # For test_socc
        T, n, m = 3, 2, 1
        self.test_qp2 = SOCP(T, n, m)
        self.z_test_socc = np.array(
            [[1.], [0.], [0.], [-3.], [-1.], [1.], [7.], [1.], [6.]])
        self.A = np.array([[1., 7.], [3., -4.]])
        self.b = np.array([[11.], [5.]])
        self.c = np.array([[3.], [7.]])
        self.d = np.array([[-1.]])
        self.AE = np.array([[.5, 7.], [3., -4.]])
        self.bE = np.array([[10.5], [5.]])
        self.cE = np.array([[3.], [6.5]])
        self.dE = np.array([[-.5]])

    def test_terms_for_socc(self):

        self.test_qp2.add_socc(self.A, self.b, self.c, self.d)
        self.test_qp2.add_socc(self.A, self.b, self.c, self.d)
        self.test_qp2.add_socc(self.AE, self.bE, self.cE, self.dE, type='end')

        x_dummy = np.array([[0], [0]])
        T, n, m = self.test_qp2.T, self.test_qp2.n, self.test_qp2.m

        self.test_qp2.P=np.zeros([0, T*(n+m)])
        P_ref = np.array([[0., 58., 128., 0., 0., 0., 0., 0., 0.],
                          [0., 0., 0., 0., 4., 212., 0., 0., 0.],
                          [0., 58., 128., 0., 0., 0., 0., 0., 0.],
                          [0., 0., 0., 0., 4., 212., 0., 0., 0.],
                          [0., 0., 0., 0., 0., 0., 0., -292., 661./2.]])

        self.assertTrue(
            (abs(self.test_qp2.P_of_zk(2*self.z_test_socc) - P_ref).sum()) < 1e-10,
            'Wrong P_of_zk(socc)') # TODO *2 weil auch im Algorithmus so, unschön

        d_ref = np.array([[-1./145.],
                          [-1./284.],
                          [-1./145.],
                          [-1./284.],
                          [-4./5371.]])

        self.assertTrue(
            (abs(self.test_qp2.form_d(x_dummy, self.z_test_socc) - d_ref).sum()) < 1e-10,
            'Wrong form_d(socc)')

        term_for_socc_ref = np.array(
            [[0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., -2./145., 52./145., 0., 0., 0., 0., 0., 0.],
             [0., 52./145., -32./145., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., -1./142., 13./71., 0., 0., 0.],
             [0., 0., 0., 0., 13./71., -8./71., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., -2./5371., 224./5371.],
             [0., 0., 0., 0., 0., 0., 0., 224./5371., -182./5371.]]) +\
                            np.array(
            [[0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., -2./145., 52./145., 0., 0., 0., 0., 0., 0.],
             [0., 52./145., -32./145., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., -1./142., 13./71., 0., 0., 0.],
             [0., 0., 0., 0., 13./71., -8./71., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0., 0., 0., 0.]])

        self.assertTrue(
            (abs(self.test_qp2.term_for_socc(self.z_test_socc) - term_for_socc_ref).sum()) < 1e-10,
            'Wrong socc_Term for Phi')

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