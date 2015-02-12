import numpy as np

import muaompc
from cvxopt import solvers, matrix as mtx
solvers.options['show_progress'] = False

def cvxopt_solve_problem(x, mpc):
    mpc.ctl.form_qp(x)
    qpx = mpc.ctl.qpx
    C = np.vstack((qpx.E, np.eye(qpx.u_lb.shape[0])))
    c_ub = np.vstack((qpx.zx_ub, qpx.u_ub))
    c_lb = np.vstack((qpx.zx_lb, qpx.u_lb))

    res = solvers.qp(mtx(qpx.HoL), mtx(qpx.gxoL),
            mtx(np.vstack((C, -C))), mtx(np.vstack((c_ub, -c_lb))))
    return np.array(res['x'])[:]

mpc = muaompc.ltidt.setup_mpc_problem('sys_motor')
x0 = np.array([-1, -1])
mpc.ctl.conf.in_iter = 24
mpc.ctl.conf.ex_iter = 2
mpc.ctl.conf.warmstart = True
xk = x0
for k in range(1):  # simulate 40 steps
    u_opt = cvxopt_solve_problem(xk, mpc)
    mpc.ctl.solve_problem(xk)
    uk = u_opt[:mpc.size.inputs, 0]
    xk = mpc.sim.predict_next_state(xk, uk)  # simulate using the optimal input

print('exact solution (cvxopt): ', u_opt)
print('approx. solution (muaompc): ', mpc.ctl.u_opt)
