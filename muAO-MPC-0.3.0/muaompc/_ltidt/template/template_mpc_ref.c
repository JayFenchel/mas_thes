#include "mpc_ref.h"
#include "mpc_inc.h"
{MPC_CONST}

/* Declaration of static functions */
void ref_compute_Rh_u_ref(struct mpc_ctl *ctl,
		real_t Rh_u_ref[], const real_t u_ref[]);
void ref_compute_Qh_x_ref(struct mpc_ctl *ctl,
		real_t Qh_x_ref[], const real_t x_ref[]);

/* Definition of external functions */

/* Compute gradient vector over Lipschitz for the given state and reference */
void ref_compute_gxoL(struct mpc_ctl *ctl, const real_t x[], 
						const real_t x_ref[], const real_t u_ref[])
{{
	{STATIC}real_t GoL_x[{HOR_INPUTS}];

	ref_compute_groL(ctl, x_ref, u_ref);
	mtx_multiply_mtx_vec(GoL_x, ctl->alm->fgm->GoL, x, {HOR_INPUTS}, {STATES});
	mtx_substract(ctl->alm->fgm->gxoL, GoL_x, ctl->alm->fgm->groL, 
					{HOR_INPUTS}, 1);
	return;
}}

/* Compute gradient vector (over Lipschitz) component for the given reference */ 
void ref_compute_groL(struct mpc_ctl *ctl,
					 const real_t x_ref[], const real_t u_ref[])
{{
	{STATIC}real_t Rh_u_ref[{HOR_INPUTS}];  /* Rh * u_ref */
	{STATIC}real_t Qh_x_ref[{HOR_STATES}];  /* Qh * x_ref */
	{STATIC}real_t Bh_T_Qh_x_ref[{HOR_INPUTS}];  /* Bh.T * Qh * x_ref */
	{STATIC}real_t gr[{HOR_INPUTS}];  /* gradient vector component of reference */

	ref_compute_Rh_u_ref(ctl, Rh_u_ref, u_ref);
	ref_compute_Qh_x_ref(ctl, Qh_x_ref, x_ref);
	mtx_multiply_mtx_vec(Bh_T_Qh_x_ref, ctl->alm->fgm->Bh_T, Qh_x_ref,
			{HOR_INPUTS}, {HOR_STATES});
	mtx_add(gr, Bh_T_Qh_x_ref, Rh_u_ref, {HOR_INPUTS}, 1);
	mtx_scale(ctl->alm->fgm->groL, gr, *(ctl->alm->Linv), {HOR_INPUTS}, 1);

	return;
}}

/* TODO make a "thread safe" version of ref_compute_groL, such that it allows
 * the computation of groL for fixed references, i.e. it will not be called
 * at every sampling time, only once a change of operation point is required.
 * _computing_ groL is time consuming, and should be carried out by a low
 * priority thread, however _setting_ it must be made _atomically_ in the 
 * MPC controller hard real-time thread. 
 */

/* Definition of static functions */

/* compute Rh * u_ref
 * Rh is a block diagonal matrix of R */
void ref_compute_Rh_u_ref(struct mpc_ctl *ctl, real_t Rh_u_ref[],
		const real_t u_ref[])
{{
	uint32_t i, ini;

	for (i=0; i<{HOR}; i++) {{
		ini = i * {INPUTS};
		mtx_multiply_mtx_vec(Rh_u_ref + ini, ctl->wmx->R, u_ref + ini,
				{INPUTS}, {INPUTS});
	}}
}}

/* compute Qh * x_ref
 * Qh is a block diagonal matrix of Q, the last block is P */
void ref_compute_Qh_x_ref(struct mpc_ctl *ctl, real_t Qh_x_ref[],
		const real_t x_ref[])
{{
	uint32_t i, ini;

	for (i=0; i<{HOR}; i++) {{
		ini = i * {STATES};
		if (i < ({HOR} - 1)) {{
			mtx_multiply_mtx_vec(Qh_x_ref + ini, ctl->wmx->Q, x_ref + ini,
					{STATES}, {STATES});
		}} else {{
			mtx_multiply_mtx_vec(Qh_x_ref + ini, ctl->wmx->P, x_ref + ini,
					{STATES}, {STATES});
		}}
	}}
}}
