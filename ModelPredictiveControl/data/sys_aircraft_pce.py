r"""
.. default-role:: math

The system description
----------------------

The system considered is the Cessna Citation 500
airfcraft presented in ([M01]_, p.64).  A continuous-time
linear model is given by `\dot{x} = A_c x + B_c u, y = C x`, where

.. math::
   A_c = \left[ \begin{matrix}
   -1.2822 & 0 & 0.98 & 0 \\
   0 & 0 & 1 & 0 \\
   -5.4293 & 0 & -1.8366 & 0 \\
   -128.2 & 128.2 & 0 & 0 \\
   \end{matrix} \right], \;\;
   B_c = \left[ \begin{matrix}
   -0.3 \\
   0 \\
   -17 \\
   0 \\
   \end{matrix} \right],
   C = \left[ \begin{matrix}
   0 & 1 & 0 & 0 \\
   0 & 0 & 0 & 1 \\
   -128.2 & 128.2 & 0 & 0 \\
   \end{matrix} \right],

and the state vector is given by `x_1` (angle of attack),
`x_2` is the pitch angle (rad), `x_3` is the pitch angle rate (rad/s), and
`x_4` is the altitude (m). The only input `u_1` is the elevator angle (rad).
The outputs are `y_1 = x_2`,  `y_2 = x_4`, and `y_3 = -128.2 x_1 + 128.2 x_2`
is the altitude rate (m/s)

The system is subject to the following constraints:

* input constraints `-0.262 \leq u_1 \leq 0.262`,
* slew rate constraint in the input `-0.524 \leq \dot{u}_1 \leq 0.524`
* state constraints `-0.349 \leq x_2 \leq 0.349`,
* output constraints `-30.0 \leq y3 \leq 30.0`.

To consider the slew rate constraint in the input, we introduce an additional
state `x_5`. The sampling interval is `dt = 0.5 s`.

The controller parameters
-------------------------

The book proposes to use identity matrices of appropriate size for
the weighting matrices `Q` and `R`. We instead select them diagonal
with values that give a similar controller performance and much lower
condition number of the Hessian of the MPC quadratic program (see
[WAD11]_), a desirable property for any numerical algorithm.

.. [M01] Maciejowski, J. M., "Predictive Control with Constraints."

.. [WAD11]  Waschl, H. and Alberer, D. and del Re, L., "Numerically Efficient Self Tuning Strategies for MPC of Integral Gas Engines".
"""

import numpy as np
from scipy import io
dt = 0.5
N = 5
mm = io.loadmat('data_matrix.mat')
Q = mm['Q_total']
R = mm['R_total']
Ad = mm['Asys']
Bd = mm['Bsys']
(n, m) = np.array(Bd).shape

P=Q
# input constraints
eui = 0.262  # rad (15 degrees). Elevator angle.
u_lb = [[-eui]]
u_ub =  [[eui]]

# mixed constraints
ex2 = 0.349  # rad/s (20 degrees). Pitch angle constraint.
ex5 = 0.524 * dt  # rad/s * dt input slew rate constraint in discrete time
ey3 = 30.
# bounds
e_lb = [[-ex2], [-ey3], [-ex5]]
e_ub = np.array([[ex2], [ex5]])
e_lb = -e_ub
(ncx, dummy) = e_ub.shape
# constraint matrices
Kx = np.zeros((ncx, n))

# terminal state constraints
data = dict(A=np.array(Ad), B=np.array(Bd), P=np.array(P), Q=np.array(Q), R=np.array(R), n=n, m=m, N=N,
u_lb=np.array(u_lb), u_ub=np.array(u_ub),
e_lb=np.array(e_lb), e_ub=np.array(e_ub),
zero=np.zeros((n,1)))
