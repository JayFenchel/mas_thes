#include "mpc_const.h" 

const real_t HoL[] = {{ /* Quadratic term (Hessian) matrix of QP over Lipschitz constant. */
{HoL} 
}};
const real_t GoL[] = {{ /* Linear term matrix of the QP, over Lipschitz constant. */
{GoL}
}};
const real_t Bh_T[] = {{ /* Extended input matrix (used for reference tracking). */
{Bh_T}
}};
const real_t E[] = {{ /* Linear factor (prediction matrix) of 2-sided state constraint. */
{E}
}};
const real_t Kx_Ai[] = {{ /* Prediction component of the state constraint bound. */
{Kx_Ai}
}};

const real_t u_lb[] = {{ /* Left (lower) constraint of the inputs for condensed QP. */
{u_lb}
}};
const real_t u_ub[] = {{ /* Right (upper) constraint of the inputs for condensed QP. */
{u_ub}
}};
const real_t e_lb[] = {{ /* Left (lower) constraint of the states for condensed QP. */
{e_lb}
}};
const real_t e_ub[] = {{ /* Right (upper) constraint of the states for condensed QP. */
{e_ub}
}};

const real_t nu = {nu}; /* Fast gradient extra step constant */
const real_t mu = {mu}; /* Augmented Lagrange multiplier penalty parameter. */
const real_t Linv = {Linv}; /* Inverse of gradient Lipschitz constant (1/L) */

/* state dependent variables */
real_t gxoL[MPC_HOR_INPUTS];  /* gradient vector as a function of the current state */
real_t zx_lb[MPC_HOR_MXCONSTRS];  /* mixed constraint lower bound as function of current state */
real_t zx_ub[MPC_HOR_MXCONSTRS];  /* mixed constraint upper bound as function of current state */

/* reference dependent variables */
real_t groL[MPC_HOR_INPUTS];

/* MPC system: state-space and weighting matrices */
const real_t Q[] = {{  /* State weighting matrix */
{Q}
}};
const real_t R[] = {{   /* Input weighting matrix */
{R}
}};
const real_t P[] = {{   /* Terminal state weighting matrix */
{P}
}};
const real_t K[] = {{   /* Linear quadratic regulator gain matrix */
{K}
}};
const real_t Ad[] = {{   /* Discrete-time system matrix */
{Ad}
}};
const real_t Bd[] = {{   /* Discrete-time input matrix */
{Bd}
}};
const real_t dt = {dt};


/* User variables created with appropriate size */
real_t u_opt[MPC_HOR_INPUTS];  /* Optimal input */
real_t l_opt[MPC_HOR_MXCONSTRS];  /* Optimal multiplier */
real_t x_trj[MPC_HOR_STATES + MPC_STATES];  /* State trajectory for u_opt */
real_t u_ini[MPC_HOR_INPUTS];  /* Initial guess input */
real_t l_ini[MPC_HOR_MXCONSTRS];  /* Initial guess multiplier */

/* Always check this declarations match the structure definitions */
struct mpc_conf conf = {{1, 1, 0}};

struct mpc_qpx qpx = {{HoL, gxoL, E, u_lb, u_ub, zx_lb, zx_ub,
		MPC_HOR_INPUTS, MPC_HOR_MXCONSTRS}};

struct mpc_sys sys = {{Ad, Bd, &dt}};

struct mpc_wmx wmx = {{Q, R, P,}};

struct mpc_lqr lqr = {{K, P,}};

struct mpc_fgm fgm = {{u_ini, gxoL, groL, &(conf.in_iter),
			HoL, GoL, Bh_T, u_lb, u_ub, &nu,
			MPC_HOR, MPC_STATES, MPC_INPUTS, MPC_HOR_INPUTS, MPC_HOR_STATES}};

struct mpc_alm alm = {{&fgm, l_ini, zx_lb, zx_ub, &(conf.ex_iter),
			&mu,
			E, Kx_Ai, e_lb, e_ub,
			&Linv, 
			MPC_STATES, MPC_MXCONSTRS, MPC_HOR_INPUTS, MPC_HOR_MXCONSTRS}};

struct mpc_ctl ctl = {{&conf, &qpx, &sys, &wmx, &lqr, &alm,
			u_opt, l_opt, x_trj, 0, 0, u_ini, l_ini}};

