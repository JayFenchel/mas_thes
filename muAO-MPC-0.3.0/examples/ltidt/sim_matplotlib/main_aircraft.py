import numpy as np

import muaompc

mpc = muaompc.ltidt.setup_mpc_problem('sys_aircraft')
# configure the controller
mpc.ctl.conf.in_iter = 24
mpc.ctl.conf.ex_iter = 2
mpc.ctl.conf.warmstart = True

# set the reference and initial value: a change of altitude
x_ref = [0, 0, 0, 400, 0] * mpc.size.hor  # a fixed reference repeated N times
mpc.ctl.x_ref = np.matrix(np.array(x_ref)).T
# mpc.ctl.u_ref is by default a zero vector
x_ini = np.zeros(mpc.size.states)
x_ini[3] = 0
mpc.sim.regulate_ref(40, x_ini)

def plot_results(data):
    plt.subplot(211)
    plt.plot(data.t, data.x[3,:])
    plt.title('Altitude [m]')
    plt.ylabel(r"$x_4$")
    plt.xlabel(r"$t$")
    plt.subplot(212)
    plt.plot(data.t, data.u[0,:], linestyle='steps-')
    plt.title('Input [rad]')
    plt.tight_layout()
    plt.show()
    plt.savefig('aircraft.eps')
    return

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
else:
    plot_results(mpc.sim.data)

