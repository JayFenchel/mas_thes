from numpy import diag  # similar effect with: from numpy import *

dt = 0.5
N = 10
mu = 100
# discrete-time system
Ad = [[  -1., 2., 0.],
      [ 0., -2., 1.],
      [ 0., 1., 1.]]
Bd = [[0, 0],
      [1, 0],
      [1, 0]]
# Weighting matrices for a problem with a better condition number
Q = diag([1., 1., 1.])
R = diag([1., 1.])
P = Q
# input constraints
eui = 0.262  # rad (15 degrees). Elevator angle.
u_lb = [[-eui],[ -eui]]
u_ub =  [[eui], [eui]]
# mixed constraints 
ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
ex5 = 0.524 * dt  # rad/s * dt input slew rate constraint in discrete time
ey3 = 30.
# bounds
e_lb = [[-ex2], [-ey3], [-ex5]]
e_ub = [[ex2], [ey3], [ex5]]
# constraint matrices
Kx = [[0, 1, 0],
      [-128.2, 128.2, 0],
      [0., 0., 0.]]
Ku = [[0, 0],
      [0, 0],
      [1, 0]]
# terminal state constraints
f_lb = e_lb
f_ub = e_ub
F = Kx
