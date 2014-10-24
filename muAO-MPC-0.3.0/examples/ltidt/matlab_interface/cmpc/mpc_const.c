#include "mpc_const.h" 

const real_t HoL[] = { /* Quadratic term (Hessian) matrix of QP over Lipschitz constant. */
0.970450065012, 0.00370019435414, 0.0035296991864, 0.00337466449212, 0.0032335773267, 0.00310506433707, 0.00298787801677, 0.00288088422044, 0.00278305081291, 0.00269343733928, 
0.00370019435414, 0.970442588812, 0.00369249730816, 0.00352175806787, 0.00336645363171, 0.0032250683554, 0.0030962259023, 0.00297867546854, 0.0028712792646, 0.00277300112789, 
0.0035296991864, 0.00369249730816, 0.970434664397, 0.00368432161003, 0.00351330465971, 0.00335769330651, 0.00321596883441, 0.00308675151199, 0.00296878678356, 0.0028609327132, 
0.00337466449212, 0.00352175806787, 0.00368432161003, 0.970426229448, 0.00367560014515, 0.00350426654549, 0.00334830524063, 0.00320619401214, 0.00307654925781, 0.00295811214406, 
0.0032335773267, 0.00336645363171, 0.00351330465971, 0.00367560014515, 0.970417211735, 0.00366625502657, 0.00349455958819, 0.0033381983897, 0.0031956452104, 0.00306551202488, 
0.00310506433707, 0.0032250683554, 0.00335769330651, 0.00350426654549, 0.00366625502657, 0.970407527324, 0.00365619563975, 0.00348408578883, 0.00332726659367, 0.00318420724976, 
0.00298787801677, 0.0030962259023, 0.00321596883441, 0.00334830524063, 0.00349455958819, 0.00365619563975, 0.970397078442, 0.00364531629967, 0.00347273071866, 0.00331538576042, 
0.00288088422044, 0.00297867546854, 0.00308675151199, 0.00320619401214, 0.0033381983897, 0.00348408578883, 0.00364531629967, 0.97038575091, 0.00363349343921, 0.00346036043548, 
0.00278305081291, 0.0028712792646, 0.00296878678356, 0.00307654925781, 0.0031956452104, 0.00332726659367, 0.00347273071866, 0.00363349343921, 0.970373411062, 0.00362058223021, 
0.00269343733928, 0.00277300112789, 0.0028609327132, 0.00295811214406, 0.00306551202488, 0.00318420724976, 0.00331538576042, 0.00346036043548, 0.00362058223021, 0.970359902037, 
 
};
const real_t GoL[] = { /* Linear term matrix of the QP, over Lipschitz constant. */
0.996336033438, 0.189731970904, 
0.993647871859, 0.180798830064, 
0.990880302337, 0.17267959956, 
0.988024973465, 0.165294954139, 
0.985072655512, 0.158572920525, 
0.982013148047, 0.152448157081, 
0.978835177851, 0.146861299853, 
0.97552628609, 0.141758368431, 
0.972072703623, 0.137090225699, 
0.968459213195, 0.132812086051, 

};
const real_t E[] = { /* Linear factor (prediction matrix) of 2-sided state constraint. */
0, 

};
const real_t Kx_Ai[] = { /* Prediction component of the state constraint bound. */
0, 

};

const real_t u_lb[] = { /* Left (lower) constraint of the inputs for condensed QP. */
-100, 
-100, 
-100, 
-100, 
-100, 
-100, 
-100, 
-100, 
-100, 
-100, 

};
const real_t u_ub[] = { /* Right (upper) constraint of the inputs for condensed QP. */
100, 
100, 
100, 
100, 
100, 
100, 
100, 
100, 
100, 
100, 

};
const real_t e_lb[] = { /* Left (lower) constraint of the states for condensed QP. */
0, 

};
const real_t e_ub[] = { /* Right (upper) constraint of the states for condensed QP. */
0, 

};

const real_t nu = 0.00847661606419; /* Fast gradient extra step constant */
const real_t mu = 0.0; /* Augmented Lagrange multiplier penalty parameter. */
const real_t Linv = 0.0; /* Inverse of gradient Lipschitz constant (1/L) */

/* state dependent variables */
real_t gxoL[MPC_HOR_INPUTS];  /* gradient vector as a function of the current state */
real_t zx_lb[MPC_HOR_STCONSTRS];  /* mixed constraint lower bound as function of current state */
real_t zx_ub[MPC_HOR_STCONSTRS];  /* mixed constraint upper bound as function of current state */

/* reference dependent variables */
real_t groL[MPC_HOR_INPUTS];

/* MPC system: state-space and weighting matrices */
const real_t Q[] = {  /* State weighting matrix */
1, 0, 
0, 1, 

};
const real_t R[] = {   /* Input weighting matrix */
1, 

};
const real_t P[] = {   /* Terminal state weighting matrix */
520.115649007, 50.0008634909, 
50.0008634909, 10.3248093315, 

};
const real_t Ad[] = {   /* Discrete-time system matrix */
1.0, 0.0095162581964, 
0.0, 0.904837418036, 

};
const real_t Bd[] = {   /* Discrete-time input matrix */
9.67483607192e-05, 
0.0190325163928, 

};
const real_t Ad_T[] = {   /* Discrete-time system matrix transpose */
1.0, 0.0, 
0.0095162581964, 0.904837418036, 

};
const real_t Bd_T[] = {   /* Discrete-time input matrix transpose */
9.67483607192e-05, 0.0190325163928, 

};
const real_t dt = 0.01;


/* User variables created with appropriate size */
real_t u_opt[MPC_HOR_INPUTS];  /* Optimal input */
real_t l_opt[MPC_HOR_STCONSTRS];  /* Optimal multiplier */
real_t u_ini[MPC_HOR_INPUTS];  /* Initial guess input */
real_t l_ini[MPC_HOR_STCONSTRS];  /* Initial guess multiplier */

/* Always check this declarations match the structure definitions */
struct mpc_conf conf = {1, 1, 0};

struct mpc_qpx qpx = {HoL, gxoL, E, u_lb, u_ub, zx_lb, zx_ub,
		MPC_HOR_INPUTS, MPC_HOR_STCONSTRS};

struct mpc_sys sys = {Ad, Bd, Ad_T, Bd_T, &dt};

struct mpc_wmx wmx = {Q, R, P,};

struct mpc_fgm fgm = {u_ini, gxoL, groL, &(conf.in_iter),
			HoL, GoL, u_lb, u_ub, &nu,
			MPC_HOR, MPC_STATES, MPC_INPUTS, MPC_HOR_INPUTS, MPC_HOR_STATES};

struct mpc_alm alm = {&fgm, l_ini, zx_lb, zx_ub, &(conf.ex_iter),
			&mu,
			E, Kx_Ai, e_lb, e_ub,
			&Linv, 
			MPC_STATES, MPC_STCONSTRS, MPC_HOR_INPUTS, MPC_HOR_STCONSTRS};

struct mpc_ctl ctl = {&conf, &qpx, &sys, &wmx, &alm,
			u_opt, l_opt, 0, 0, u_ini, l_ini};

