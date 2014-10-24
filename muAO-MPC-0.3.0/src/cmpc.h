#ifndef MPCPY_CMPC_H
#define MPCPY_CMPC_H

#include <Python.h>
#include "mpc.h"
#include <stdio.h>

void python_add_matrix(
	PyObject* dict,
	const char* name,
	const real_t* data,
	uint32_t rows,
	uint32_t cols
);

void python_read_matrix(
	PyObject* structure,
	const char* name,
	real_t* data,
	uint32_t rows,
	uint32_t cols
);

PyObject* python_create_conf(const struct mpc_ctl* data);
PyObject* python_create_alm(const struct mpc_ctl* data);
PyObject* python_create_fgm(const struct mpc_ctl* data);
PyObject* python_create_qpx(const struct mpc_ctl* data);
PyObject* python_create_wmx(const struct mpc_ctl* data);
PyObject* python_create_lqr(const struct mpc_ctl* data);
PyObject* python_create_sys(const struct mpc_ctl* data);
PyObject* python_create_ctl(const struct mpc_ctl* data);

#endif
