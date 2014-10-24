#ifndef MPC_CONST_H
#define MPC_CONST_H

#include "mpc_base.h"



enum { 
MPC_HOR = 10,  /**< MPC prediction horizon. */
MPC_STATES = 2,  /**< Number of system states. */
MPC_INPUTS = 1,  /**< Number of system inputs. */
MPC_STCONSTRS = 0, /**< Number of mixed stage constraints. */
MPC_HOR_INPUTS = 10,  /**< Horizon times number of inputs. */
MPC_HOR_STATES = 20,  /**< Horizon times number of states. */
MPC_HOR_STCONSTRS = 0,  /**< Horizon times number of mixed constrained
 	 	 	 	 	 	 	 	 	plus the number of end state constraints. */
}; 

#endif /* MPC_CONST_H */
