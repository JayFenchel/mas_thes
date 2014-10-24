#!/usr/bin/python3

import numpy as np
from scipy.io import savemat
import matplotlib as mpl
import matplotlib.pyplot as plt

from cvxopt import solvers, matrix as mtx
solvers.options['show_progress'] = False

import muaompc
x0 = np.array([0., 0., 0., 50., 0.])  # simulation initial state
points = 20  # number of simulation steps
fac = 1e2  # factor for initial guesses

def main():
    import sys_aircraft as sys
    n, m = np.array(sys.Bd).shape
    # Original weighting matrices
    sys.Q = np.eye(n)
    sys.R = np.eye(m)
    sys.P = sys.Q
    mpc = muaompc.ltidt.setup_mpc_problem(sys)
    print('Original problem:')
    print('Internal c.n.: ', mpc.tun.calc_int_cn()[1])
    print('Hessian c.n.: ', np.linalg.cond(mpc.ctl.qpx.HoL))
    data = simulate_controller(mpc)
    Qd, Rd = find_wmx_diag(mpc, data)
    print('Qd:', Qd)
    print('Rd:', Rd)
    # weighting matrices for lower condition number
    sys.Q = np.diag(Qd)
    sys.R = np.diag(Rd)
    sys.P = sys.Q
    mpc = muaompc.ltidt.setup_mpc_problem(sys)
    print('Conditioned problem:')
    print('Internal c.n.: ', mpc.tun.calc_int_cn()[1])
    print('Hessian c.n.: ', np.linalg.cond(mpc.ctl.qpx.HoL))
    write_wmx(sys.Q, sys.R)

def simulate_controller(mpc):
    def solve_problem_cvxopt(x, mpc=mpc):
        mpc.ctl.form_qp(x)
        qpx = mpc.ctl.qpx
        C = np.vstack((qpx.E, np.eye(qpx.u_lb.shape[0])))
        c_ub = np.vstack((qpx.zx_ub, qpx.u_ub))
        c_lb = np.vstack((qpx.zx_lb, qpx.u_lb))

        res = solvers.qp(mtx(qpx.HoL), mtx(qpx.gxoL),
                mtx(np.vstack((C, -C))), mtx(np.vstack((c_ub, -c_lb))))
        mpc.ctl.u_opt[:] = np.array(res['x'])[:]
        return
    mpc.ctl.solve_problem = solve_problem_cvxopt
    mpc.sim.regulate_ref(points, x0)
    return mpc.sim.data

def find_wmx_diag(mpc, data):
    xref = data.x
    uref = data.u0
    stw_mpc = np.ones(xref.shape[0]) * fac
    inw_mpc = np.ones(uref.shape[0]) * fac
    r = mpc.tun.reduce_H_cn(xref, uref, np.diag(mpc.wmx.Q),
            np.diag(mpc.wmx.R), cn_ub=25, stw0=stw_mpc, inw0=inw_mpc)
    Q = r['stw']
    R = r['inw']
    return Q, R

def write_wmx(Q, R):
    with open('wmx.py', 'w') as f:
        f.write('Q = ' + str(Q.tolist()) + '\n')
        f.write('R = ' + str(R.tolist()) + '\n')

if __name__ == '__main__':
    main()
