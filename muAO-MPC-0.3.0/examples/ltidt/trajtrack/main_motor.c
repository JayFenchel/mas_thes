#include <stdio.h>  /* printf */

#include "cmpc/include/mpc.h"  /* the auto-generated code */

/* This file is a test of the C routines of the ALM+FGM MPC
 * algorithm. The same routines can be used in real systems.
 * It demonstrates the syntax of reference tracking problems.
 */
int main(void)
{
	extern struct mpc_ctl ctl;  /* already defined */
	real_t x[MPC_STATES];  /* current state of the system */
	real_t u_ref[MPC_HOR_INPUTS];  /* the input reference */
	real_t x_ref[MPC_HOR_STATES];  /* the state reference */

	int i;

	/* set the input and state reference to the desired value */	
	for (i=0; i<MPC_HOR_INPUTS; i++) {
		u_ref[i] = 0.; 
	}
	for (i=0; i<MPC_HOR_STATES; i++) {
		x_ref[i] = 1.;
	}

	/* assign the references to the mpc_ctl respective pointers */
	ctl.u_ref = u_ref;
	ctl.x_ref = x_ref;

	ctl.conf->in_iter = 15;  /* number of iterations */

	/* The current state */
	x[0] = 0.;
	x[1] = 0.;

	/* Solve MPC problem and print the input sequence */
	mpc_ctl_solve_problem(&ctl, x); 
	for (i=0; i<MPC_HOR_INPUTS; i++) {
		printf("u[0] = %f \n", (ctl.u_opt[i]));
	}

	return 0;
}
