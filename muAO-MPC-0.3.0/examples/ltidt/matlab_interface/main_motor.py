"""Create the MPC controller MATLAB interface.

For the complete example, run the main_motor.m MATLAB script.
"""

import muaompc

mpc = muaompc.ltidt.setup_mpc_problem('sys_motor.mat')
mpc.generate_c_files(matlab=True)

