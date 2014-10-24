r"""Example use of the trajectory tracking problem.

By default, the MPC controller stabilizes system at the origin.
If the controller is to stabilize a different point, or 
track a trajectory is defined online, not at this stage.

See the file main_motor.c for an example of how to define the 
trajectory to track.
"""
import muaompc
mpc = muaompc.ltidt.setup_mpc_problem('sys_legomotor')
mpc.generate_c_files()
