/** Solve state constrained (stc) Model Predictive Control (MPC) problem.
 *
 */

#include <string.h> /* sizeof */
#include "mpc_inc.h"
#include "mpc_stc.h"
{MPC_CONST}

/* Declaration of static functions */

static void stc_fgm_minimize_qp(const struct mpc_alm *alm, real_t u[],
		const real_t u_i[], const real_t l_i[]);
static void stc_fgm_compute_grad_over_L(const struct mpc_alm *alm,
		real_t gradoL[], const real_t w[], const real_t l_i[]);
static void stc_fgm_compute_state_constr_violation_penalty(
		const struct mpc_alm *alm, real_t out[], const real_t w[],
		const real_t l_i[], const real_t zx_b[]);
static void stc_fgm_compute_l_update(const struct mpc_alm *alm,
		real_t l_i1[], const real_t w[], const real_t l_i[]);
static void stc_fgm_compute_grad_alm(const struct mpc_alm *alm,
		real_t out[], const real_t w[],	const real_t l_i[]);


/* Definition of external function */

/* Form quadratic program for current state */
void stc_ctl_form_qp(struct mpc_ctl *ctl, const real_t x[])
{{
	inc_ctl_form_qp(ctl, x);
	stc_compute_state_constr_bound_online(ctl->alm, ctl->alm->zx_lb, x, 
			ctl->alm->e_lb);
	stc_compute_state_constr_bound_online(ctl->alm, ctl->alm->zx_ub, x, 
			ctl->alm->e_ub);
	
	return;
}}

/* Warmstart state constrained MPC */
void stc_ctl_warmstart(struct mpc_ctl *ctl) {{
	inc_ctl_warmstart(ctl);
	mpc_warmstart_vector(ctl->alm->fgm, ctl->alm->l_0, ctl->l_opt, {MXCONSTRS});
	
	return;
}}

/* Solve state constrained MPC problem for the current state */
void stc_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[])
{{
	stc_ctl_form_qp(ctl, x);
	stc_alm_minimize_qp(ctl->alm, ctl->u_opt, ctl->l_opt);
	
	return;
}}

/* Minimize MPC quadratic program using augmented Lagrangian method
 * together with fast gradient method */
void stc_alm_minimize_qp(struct mpc_alm *alm,
		real_t u_opt[], real_t l_opt[])
{{

	{STATIC}real_t u_i[{HOR_INPUTS}]; /* iteration update of input */
	{STATIC}real_t l_i[{HOR_MXCONSTRS}]; /* it. update of Lagrange mult.*/

	uint32_t i;
	const size_t sizeof_u = {HOR_INPUTS} * sizeof(*u_opt);
	const size_t sizeof_l = {HOR_MXCONSTRS} * sizeof(*l_opt);

	memcpy(u_i, alm->fgm->u_0, sizeof_u);
	memcpy(l_i, alm->l_0, sizeof_l);

	for (i = 0; i < *(alm->i_ex); i++) {{
		stc_fgm_minimize_qp(alm, u_opt, u_i, l_i);
		stc_fgm_compute_l_update(alm, l_opt, u_opt, l_i);
		/* do not take the next two outside the loop. The two previous
		 * functions expect different pointers for outputs u/l_upd
		 * and inputs u_i/l_i respectively. */
		memcpy(u_i, u_opt, sizeof_u);
		memcpy(l_i, l_opt, sizeof_l);
	}}

	return;
}}


void stc_compute_state_constr_bound_online(const struct mpc_alm *alm,
		real_t out[], const real_t x[], const real_t bound[])
{{
	{STATIC}real_t scb_x[{HOR_MXCONSTRS}];

	mtx_multiply_mtx_vec(scb_x, alm->Kx_Ai, x, {HOR_MXCONSTRS}, {STATES});
	mtx_substract(out, bound, scb_x, {HOR_MXCONSTRS}, 1);

	return;
}}


/* Definition of static functions */

/* Minimize ALM internal iteration using a fast gradient method */
static void stc_fgm_minimize_qp(const struct mpc_alm *alm, real_t u[],
		const real_t u_i[], const real_t l_i[])
{{
	{STATIC}real_t u_old[{HOR_INPUTS}];
	{STATIC}real_t w[{HOR_INPUTS}];
	{STATIC}real_t gradoL[{HOR_INPUTS}];
	uint32_t i;
	const size_t sizeof_u = {HOR_INPUTS} * sizeof(*u);

	memcpy(u, u_i, sizeof_u);
	memcpy(w, u, sizeof_u);
	memcpy(u_old, u, sizeof_u);

	for (i = 0; i < *(alm->fgm->j_in); i++) {{
		stc_fgm_compute_grad_over_L(alm, gradoL, w, l_i);
		inc_fgm_minimize_qp_iteration(alm->fgm, u, u_old, w, gradoL);
	}}

	return;
}}

/* Compute gradient divided by Lipschitz constant */
static void stc_fgm_compute_grad_over_L(const struct mpc_alm *alm,
		real_t gradoL[], const real_t w[], const real_t l_i[])
{{
	{STATIC}real_t gradoL_inc[{HOR_INPUTS}]; /* input constrained grad. over L*/
	{STATIC}real_t grad_stc[{HOR_INPUTS}]; /* state constrained grad. */
	{STATIC}real_t gradoL_stc[{HOR_INPUTS}]; /* state constrained grad. over L*/

	inc_fgm_compute_grad_over_L(alm->fgm, gradoL_inc, w);
	stc_fgm_compute_grad_alm(alm, grad_stc, w, l_i);
	mtx_scale(gradoL_stc, grad_stc, *(alm->Linv), {HOR_INPUTS}, 1);
	mtx_add(gradoL, gradoL_inc, gradoL_stc, {HOR_INPUTS}, 1);

	return;
}}

static void stc_fgm_compute_state_constr_violation_penalty(
		const struct mpc_alm *alm, real_t out[], const real_t w[],
		const real_t l_i[], const real_t zx_b[])
{{
	/* penalty = l_i + mu * (E * w - z_b) */
	{STATIC}real_t E_w[{HOR_MXCONSTRS}];
	{STATIC}real_t diff[{HOR_MXCONSTRS}];
	{STATIC}real_t pen_diff[{HOR_MXCONSTRS}];

	mtx_multiply_mtx_vec(E_w, alm->E, w, {HOR_MXCONSTRS}, {HOR_INPUTS});
	mtx_substract(diff, E_w, zx_b, {HOR_MXCONSTRS}, 1);
	mtx_scale(pen_diff, diff, *(alm->mu), {HOR_MXCONSTRS}, 1);
	mtx_add(out, l_i, pen_diff, {HOR_MXCONSTRS}, 1);

	return;
}}

/* Compute multiplier update */
static void stc_fgm_compute_l_update(const struct mpc_alm *alm,
		real_t l_i1[], const real_t w[], const real_t l_i[])
{{
	{STATIC}real_t p_pos[{HOR_MXCONSTRS}]; /* positive penalty */
	{STATIC}real_t p_neg[{HOR_MXCONSTRS}]; /* negative penalty */

	stc_fgm_compute_state_constr_violation_penalty(alm, p_pos, w, l_i, alm->zx_ub);
	stc_fgm_compute_state_constr_violation_penalty(alm, p_neg, w, l_i, alm->zx_lb);

	mtx_max_vec_zero(p_pos, {HOR_MXCONSTRS});
	mtx_min_vec_zero(p_neg, {HOR_MXCONSTRS});

	mtx_add(l_i1, p_pos, p_neg, {HOR_MXCONSTRS}, 1);

	return;
}}

static void stc_fgm_compute_grad_alm(const struct mpc_alm *alm,
		real_t out[], const real_t w[],	const real_t l_i[])
{{
	{STATIC}real_t l_grad[{HOR_MXCONSTRS}];
	{STATIC}real_t E_T[{HOR_MXCONSTRS} * {HOR_INPUTS}];

	stc_fgm_compute_l_update(alm, l_grad, w, l_i);
	mtx_transpose(E_T, alm->E, {HOR_MXCONSTRS}, {HOR_INPUTS});
	mtx_multiply_mtx_vec(out, E_T, l_grad, {HOR_INPUTS}, {HOR_MXCONSTRS});

	return;
}}
