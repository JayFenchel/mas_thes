extern "C" {
#include "cmpc/include/mpc.h"
}
#include <qpOASES.hpp>

/* This example demonstrates the use of muAO-MPC with a different QP solver.
 * In this case, we use qpOASES-3.0Beta, a powerful QP solver for MPC.
 */
int main(void)
{
	extern struct mpc_ctl ctl;  /* the problem data */
	struct mpc_qpx *qpx = ctl.qpx;  /* a pointer copy, just for convenience */

	/* the qpOASES initialization */
	USING_NAMESPACE_QPOASES
	QProblem qp(MPC_HOR_INPUTS, MPC_HOR_MXCONSTRS);
	Options opt;
	opt.setToFast();
	opt.printLevel = PL_LOW;
	qp.setOptions(opt);
	int nWSR = 15;
	returnValue ret;
	real_t up[MPC_HOR_INPUTS];  /* primal solution */

	real_t x[MPC_STATES] = {0., 0., 0., -400., 0.};  /* initial state */
	mpc_ctl_form_qp(&ctl, x);  /* form a QP for x, and store it in ctl.qpx */
	/* solve for initial state using qpOASES */
	ret = qp.init(qpx->HoL, qpx->gxoL, 
			qpx->E, qpx->u_lb, qpx->u_ub,
			qpx->zx_lb, qpx->zx_ub,
			nWSR, 0);


	if (SUCCESSFUL_RETURN == ret) {
		qp.getPrimalSolution(up);

		for (int j=0; j<MPC_HOR_INPUTS; j++) {
			printf("u[%d]=%e ", j, up[j]);
		}
	}
	printf("\n");
	/* the new state comes from simulation in this case */
	x[0]=2.469694e-01;
	x[1]=2.877036e-01;
	x[2]=8.967119e-01;
	x[3]=-3.996400e+02;
	x[4]=-2.000328e-01;

	/* solve again for new state using a hot-start */
	mpc_ctl_form_qp(&ctl, x);
	ret = qp.hotstart(qpx->gxoL, qpx->u_lb, qpx->u_ub,
			qpx->zx_lb, qpx->zx_ub, nWSR, 0);

	if (SUCCESSFUL_RETURN == ret) {
		qp.getPrimalSolution(up);

		for (int j=0; j<MPC_HOR_INPUTS; j++) {
			printf("u[%d]=%e ", j, up[j]);
		}
	}
	printf("\n");

	return 0;
}
