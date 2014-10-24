#include <stdio.h>  /* printf */

#include "cmpc/include/mpc.h"  /* the auto-generated code */
#ifdef FIP_OPS  /* defined automatically by the code generation */
#include "cmpc/include/fip_ops.h"
#endif
/* This file is a test of the C routines of the ALM+FGM MPC
 * algorithm. The same routines can be used in real systems.
 * It uses a experimental fixed-point arithmetic module.
 */
int main(void)
{
	real_t x[MPC_STATES];  /* current state of the system */
	const float32_t X0 = 2., X1 = 14.;  /* sample state */
	extern struct mpc_ctl ctl;  /* already defined */

	ctl.conf->in_iter = 10;  /* number of iterations */

	/* The current state */
#ifdef FIP_OPS  /* transform a 32-bit float into 32-bit fixed-point */
	x[0] = fip_real2fip(X0);
	x[1] = fip_real2fip(X1);
#else
	x[0] = X0;
	x[1] = X1;
#endif
	/* Solve MPC problem and print the first element of input sequence */
	mpc_ctl_solve_problem(&ctl, x);  /* solve the MPC problem */
#ifdef FIP_OPS
	printf("u[0] = %f \n", fip_fip2real(ctl.u_opt[0]));
#else
	printf("u[0] = %f \n", (ctl.u_opt[0]));
#endif
	printf("Compare to the double precision solution computed in Python: \n");
	printf("u[0] = -35.13150044 \n");

	return 0;
}
