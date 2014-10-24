#ifndef MPC_CONST_H
#define MPC_CONST_H

#include "mpc_base.h"

{extra_defines}

enum {{ 
MPC_HOR = {HOR},  /**< MPC prediction horizon. */
MPC_STATES = {STATES},  /**< Number of system states. */
MPC_INPUTS = {INPUTS},  /**< Number of system inputs. */
MPC_MXCONSTRS = {MXCONSTRS}, /**< Number of mixed stage constraints. */
MPC_HOR_INPUTS = {HOR_INPUTS},  /**< Horizon times number of inputs. */
MPC_HOR_STATES = {HOR_STATES},  /**< Horizon times number of states. */
MPC_HOR_MXCONSTRS = {HOR_MXCONSTRS}  /**< Horizon times number of mixed constrained
plus the number of end state constraints. */
}}; 

#endif /* MPC_CONST_H */
