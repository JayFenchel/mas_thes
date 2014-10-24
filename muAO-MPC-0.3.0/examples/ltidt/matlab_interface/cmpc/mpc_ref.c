#include "mpc_ref.h"
#include "mpc_inc.h"
#include "mpc_const.h"

/* Declaration of static functions */
static void ref_compute_ref_qp_grad_vec(struct mpc_ctl *ctl, real_t Fref[],
		const real_t stwmref[], const real_t inwmref[], const real_t tewmref[]);
static void ref_compute_endref_weight_mtx(struct mpc_ctl *ctl, real_t wmtxref[],
		const real_t vecref[], const real_t wmtxbase[], const int vec_basesize);
static void ref_compute_ref_weight_mtx(struct mpc_ctl *ctl, real_t wmtxref[], 
    const real_t vecref[], const real_t wmtxbase[], const int vec_basesize);
static void ref_compute_endref_qp_grad_vec(struct mpc_ctl *ctl,
		real_t AdT_pow_k[], real_t Fref_end[],
		const real_t stwmref[], const real_t inwmref[], const real_t tewmref[]);

/* Definition of external functions */

/* Compute gradient vector over Lipschitz for the given state and reference */
void ref_compute_gxoL(struct mpc_ctl *ctl, const real_t x[], 
						const real_t x_ref[], const real_t u_ref[])
{
	static real_t GoL_x[MPC_HOR_INPUTS];

	ref_compute_groL(ctl, x_ref, u_ref);
	mtx_multiply_mtx_vec(GoL_x, ctl->alm->fgm->GoL, x, MPC_HOR_INPUTS, MPC_STATES);
	mtx_substract(ctl->alm->fgm->gxoL, GoL_x, ctl->alm->fgm->groL, 
					MPC_HOR_INPUTS, 1);

	return;
}

/* Compute gradient vector (over Lipschitz) component for the given reference */ 
void ref_compute_groL(struct mpc_ctl *ctl,
					 const real_t x_ref[], const real_t u_ref[])
{
	static real_t gr[MPC_HOR_INPUTS];  /* gradient vector component of reference */
	static real_t stwmref[MPC_HOR_STATES];  /* state reference weight matrix */
	static real_t inwmref[MPC_HOR_INPUTS];  /* input reference weight matrix */
	static real_t tewmref[MPC_STATES];  /* terminal state reference weight matrix */

	ref_compute_ref_weight_mtx(ctl, stwmref, x_ref, ctl->wmx->Q, MPC_STATES);
	ref_compute_ref_weight_mtx(ctl, inwmref, u_ref, ctl->wmx->R, MPC_INPUTS);
	ref_compute_endref_weight_mtx(ctl, tewmref, x_ref, ctl->wmx->P, MPC_STATES);

	ref_compute_ref_qp_grad_vec(ctl, gr, stwmref, inwmref, tewmref);
	mtx_scale(ctl->alm->fgm->groL, gr, *(ctl->alm->Linv), MPC_HOR_INPUTS, 1);

	return;
}

/* TODO make a "thread safe" version of ref_compute_groL, such that it allows
 * the computation of groL for fixed references, i.e. it will not be called
 * at every sampling time, only once a change of operation point is required.
 * _computing_ groL is time consuming, and should be carried out by a low
 * priority thread, however _setting_ it must be made _atomically_ in the 
 * MPC controller hard real-time thread. 
 */

/* Definition of static functions */

static void ref_compute_endref_weight_mtx(struct mpc_ctl *ctl, real_t wmtxref[],
							const real_t vecref[],
							const real_t wmtxbase[], const int vec_basesize)
{
    const real_t *pvecend = vecref + (MPC_HOR - 1) * vec_basesize;
    mtx_multiply_mtx_vec(wmtxref, wmtxbase, pvecend,
    						vec_basesize, vec_basesize);

    return;
}

static void ref_compute_ref_weight_mtx(struct mpc_ctl *ctl, real_t wmtxref[],
							const real_t vecref[],
							const real_t wmtxbase[], const int vec_basesize)
{
  unsigned i, ini;
  real_t *pwmtxref_i;
  const real_t *pvecref_i;

  for (i=0; i < MPC_HOR; i++) {
    ini = i * vec_basesize;
    pwmtxref_i = wmtxref + ini;
    pvecref_i = vecref + ini;
    mtx_multiply_mtx_vec(pwmtxref_i, wmtxbase, pvecref_i,
    					vec_basesize, vec_basesize);
  }
  return;
}

static void ref_compute_endref_qp_grad_vec(struct mpc_ctl *ctl,
					real_t AdT_pow_k[], real_t Fref_end[],
					const real_t stwmref[], const real_t inwmref[],
					const real_t tewmref[])
{
  struct mpc_sys *sys = ctl->sys;
	unsigned in_ini, st_ini, i;
	real_t Ad_T_tewmref[MPC_STATES];
	real_t Bd_T_tewmref[MPC_INPUTS];
	const real_t *pstwmref_end, *pinwmref_end;

	i = MPC_HOR - 1;
	in_ini = i * MPC_INPUTS;
        st_ini = i * MPC_STATES;

	pstwmref_end = stwmref + st_ini;
	mtx_multiply_mtx_vec(Ad_T_tewmref, sys->Ad_T, tewmref, MPC_STATES, MPC_STATES);
	mtx_add(AdT_pow_k, Ad_T_tewmref, pstwmref_end, MPC_STATES, 1);

	pinwmref_end = inwmref + in_ini;
	mtx_multiply_mtx_vec(Bd_T_tewmref, sys->Bd_T, tewmref, MPC_INPUTS, MPC_STATES);
	mtx_add(Fref_end, Bd_T_tewmref, pinwmref_end, MPC_INPUTS, 1);

	return;
}


static void ref_compute_ref_qp_grad_vec(struct mpc_ctl *ctl,
		real_t Fref[],
					const real_t stwmref[],
					const real_t inwmref[],
					const real_t tewmref[])
{
  struct mpc_sys *sys = ctl->sys;
  int in_ini, st_ini, i;

  static real_t AdT_pow_k[MPC_STATES];
  static real_t BdT_AdTk[MPC_INPUTS];
  static real_t AdT_AdTk[MPC_STATES];
  real_t *pFref_k;
  real_t *pFref_end = Fref + (MPC_HOR - 1) * MPC_INPUTS;

  ref_compute_endref_qp_grad_vec(ctl, AdT_pow_k, pFref_end, 
		  stwmref, inwmref, tewmref);

  for (i=(MPC_HOR - 2); i > -1; i--) {

     mtx_multiply_mtx_vec(BdT_AdTk, sys->Bd_T, AdT_pow_k, MPC_INPUTS, MPC_STATES);
     in_ini = i * MPC_INPUTS;
     pFref_k = Fref + in_ini;

     mtx_add(pFref_k, BdT_AdTk, inwmref + in_ini, MPC_INPUTS, 1);
     mtx_multiply_mtx_vec(AdT_AdTk, sys->Ad_T, AdT_pow_k, MPC_STATES, MPC_STATES);
     st_ini = i * MPC_STATES;

     mtx_add(AdT_pow_k, AdT_AdTk, stwmref + st_ini, MPC_STATES, 1);
  }
}
