#include <stddef.h>  /* NULL */
#include "mpc.h"
#include "mpc_const.h"

static void inc_ref_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[]);
static void stc_ref_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[]);

void mpc_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[])
{
	if ((ctl->x_ref == NULL) || (ctl->u_ref == NULL)) { /* any pointer is zero */
		if (MPC_HOR_STCONSTRS == 0) {
			inc_ctl_solve_problem(ctl, x);
		} else {
			stc_ctl_solve_problem(ctl, x);
		}
	} else { /* both references are given */
		if (MPC_HOR_STCONSTRS == 0) {
			inc_ref_ctl_solve_problem(ctl, x);
		} else {
			stc_ref_ctl_solve_problem(ctl, x);
		}
	}

	if (ctl->conf->warmstart) {
		if (MPC_HOR_STCONSTRS == 0) {
			inc_ctl_warmstart(ctl);
		} else {
			stc_ctl_warmstart(ctl);
		}
	}
	
	return;
}

void mpc_ctl_form_qp(struct mpc_ctl *ctl, const real_t x[])
{
	const uint32_t in_iter = ctl->conf->in_iter;
	const uint32_t ex_iter = ctl->conf->ex_iter;

	/* Form the appropriate QP without actually solving the MPC problem */
	ctl->conf->in_iter = 0;
	ctl->conf->ex_iter = 0;
	mpc_ctl_solve_problem(ctl, x);
	
	ctl->conf->in_iter = in_iter;
	ctl->conf->ex_iter = ex_iter;
	
	return;
}

void mpc_predict_next_state(struct mpc_ctl *ctl, real_t x[])
{
	real_t Adx[MPC_STATES], Bdu[MPC_STATES];

	/* compute x = Ad * x + Bd * u_opt; it overwrites x */
	mtx_multiply_mtx_vec(Adx, ctl->sys->Ad, x, MPC_STATES, MPC_STATES);
	mtx_multiply_mtx_vec(Bdu, ctl->sys->Bd, ctl->u_opt, MPC_STATES, MPC_INPUTS);
	mtx_add(x, Adx, Bdu, MPC_STATES, 1);

	return;
}

static void inc_ref_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[])
{
	ref_compute_gxoL(ctl, x, ctl->x_ref, ctl->u_ref);
	inc_fgm_minimize_qp(ctl->alm->fgm, ctl->u_opt);
	
	return;
}

static void stc_ref_ctl_solve_problem(struct mpc_ctl *ctl, const real_t x[])
{
	ref_compute_gxoL(ctl, x, ctl->x_ref, ctl->u_ref);
	stc_alm_minimize_qp(ctl->alm, ctl->u_opt, ctl->l_opt);
	
	return;
}
