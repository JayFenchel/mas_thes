#include "cmpc.h"

#define PY_ARRAY_UNIQUE_SYMBOL _muaompc_cmpc_numpy
#define NO_IMPORT_ARRAY
#include <numpy/npy_common.h>
#include <numpy/arrayobject.h>

void python_add_matrix(
	PyObject* dict,
	const char* name,
	const real_t* data,
	uint32_t rows,
	uint32_t cols
) {
	npy_intp dims[2];
	dims[0] = rows;
	dims[1] = cols;

	PyArrayObject* matrix = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_FLOAT64);
	memcpy(PyArray_DATA(matrix), data, sizeof(float64_t) * dims[0] * dims[1]);

	PyDict_SetItemString(dict, name, (PyObject*)matrix);
}

void python_read_matrix(
	PyObject* structure,
	const char* name,
	real_t* data,
	uint32_t rows,
	uint32_t cols
) {
	PyArrayObject* matrix = (PyArrayObject*)PyDict_GetItemString(structure, name);
	memcpy(data, PyArray_DATA(matrix), sizeof(float64_t) * rows * cols);
}

PyObject* python_create_conf(const struct mpc_ctl* data) {
	PyObject* conf = PyDict_New();

	PyObject* in_iter = PyLong_FromLong(data->conf->in_iter);
	PyObject* ex_iter = PyLong_FromLong(data->conf->ex_iter);
	PyObject* warmstart = PyBool_FromLong(data->conf->warmstart);

	PyDict_SetItemString(conf, "ex_iter", ex_iter);
	PyDict_SetItemString(conf, "in_iter", in_iter);
	PyDict_SetItemString(conf, "warmstart", warmstart);

	return conf;
}

PyObject* python_create_sys(const struct mpc_ctl* data) {
	PyObject* sys = PyDict_New();

	python_add_matrix(
		sys, "Ad", data->sys->Ad,
		data->alm->STATES, data->alm->STATES
	);

	python_add_matrix(
		sys, "Bd", data->sys->Bd,
		data->alm->fgm->STATES, data->alm->fgm->INPUTS
	);

	PyObject* dt = PyFloat_FromDouble(*(data->sys->dt));
	PyDict_SetItemString(sys, "dt", dt);

	return sys;
}

PyObject* python_create_wmx(const struct mpc_ctl *data) {
	PyObject* wmx = PyDict_New();

	python_add_matrix(
		wmx, "Q", data->wmx->Q,
		data->alm->fgm->STATES, data->alm->fgm->STATES
	);

	python_add_matrix(
		wmx, "R", data->wmx->R,
		data->alm->fgm->INPUTS, data->alm->fgm->INPUTS
	);

	python_add_matrix(
		wmx, "P", data->wmx->P,
		data->alm->fgm->STATES, data->alm->fgm->STATES
	);

	return wmx;
}

PyObject* python_create_lqr(const struct mpc_ctl *data) {
	PyObject* lqr = PyDict_New();

	python_add_matrix(
		lqr, "K", data->lqr->K,
		data->alm->fgm->INPUTS, data->alm->fgm->STATES
	);

	python_add_matrix(
		lqr, "P", data->lqr->P,
		data->alm->fgm->STATES, data->alm->fgm->STATES
	);

	return lqr;
}

PyObject* python_create_qpx(const struct mpc_ctl* data) {
	PyObject* qpx = PyDict_New();

	python_add_matrix(
		qpx, "HoL", data->qpx->HoL,
		data->qpx->HOR_INPUTS, data->qpx->HOR_INPUTS
	);

	python_add_matrix(
		qpx, "gxoL", data->qpx->gxoL,
		data->qpx->HOR_INPUTS, 1
	);

	python_add_matrix(
		qpx, "E", data->qpx->E,
		data->qpx->HOR_MXCONSTRS, data->qpx->HOR_INPUTS
	);

	python_add_matrix(
		qpx, "u_lb", data->qpx->u_lb,
		data->qpx->HOR_INPUTS, 1
	);

	python_add_matrix(
		qpx, "u_ub", data->qpx->u_ub,
		data->qpx->HOR_INPUTS, 1
	);

	python_add_matrix(
		qpx, "zx_lb", data->qpx->zx_lb,
		data->qpx->HOR_MXCONSTRS, 1
	);

	python_add_matrix(
		qpx, "zx_ub", data->qpx->zx_ub,
		data->qpx->HOR_MXCONSTRS, 1
	);

	PyObject* HOR_INPUTS = PyLong_FromLong(data->qpx->HOR_INPUTS);
	PyDict_SetItemString(qpx, "HOR_INPUTS", HOR_INPUTS);

	PyObject* HOR_MXCONSTRS = PyLong_FromLong(data->qpx->HOR_MXCONSTRS);
	PyDict_SetItemString(qpx, "HOR_MXCONSTRS", HOR_MXCONSTRS);

	return qpx;
}

PyObject* python_create_fgm(const struct mpc_ctl* data) {
	PyObject* fgm = PyDict_New();

	python_add_matrix(
		fgm, "u_0", data->alm->fgm->u_0,
		data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(
		fgm, "gxoL", data->alm->fgm->gxoL,
		data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(
		fgm, "groL", data->alm->fgm->groL,
		data->alm->fgm->HOR_INPUTS, 1
	);

	PyObject* j_in = PyLong_FromLong(*(data->alm->fgm->j_in));
	PyDict_SetItemString(fgm, "j_in", j_in);

	python_add_matrix(
		fgm, "HoL", data->alm->fgm->HoL,
		data->alm->fgm->HOR_INPUTS,
		data->alm->fgm->HOR_INPUTS
	);

	python_add_matrix(
		fgm, "GoL", data->alm->fgm->GoL,
		data->alm->fgm->HOR_INPUTS,
		data->alm->fgm->STATES
	);

	python_add_matrix(
		fgm, "Bh_T", data->alm->fgm->Bh_T,
		data->alm->fgm->HOR_INPUTS,
		data->alm->fgm->HOR_STATES
	);

	python_add_matrix(
		fgm, "u_lb", data->alm->fgm->u_lb,
		data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(
		fgm, "u_ub", data->alm->fgm->u_ub,
		data->alm->fgm->HOR_INPUTS, 1
	);

	PyObject* nu = PyFloat_FromDouble(*(data->alm->fgm->nu));
	PyDict_SetItemString(fgm, "nu", nu);

	PyObject* STATES = PyLong_FromLong(data->alm->fgm->STATES);
	PyDict_SetItemString(fgm, "STATES", STATES);

	PyObject* INPUTS = PyLong_FromLong(data->alm->fgm->INPUTS);
	PyDict_SetItemString(fgm, "INPUTS", INPUTS);

	PyObject* HOR_INPUTS = PyLong_FromLong(data->alm->fgm->HOR_INPUTS);
	PyDict_SetItemString(fgm, "HOR_INPUTS", HOR_INPUTS);

	PyObject* HOR = PyLong_FromLong(data->alm->fgm->HOR);
	PyDict_SetItemString(fgm, "HOR", HOR);

	PyObject* HOR_STATES = PyLong_FromLong(data->alm->fgm->HOR_STATES);
	PyDict_SetItemString(fgm, "HOR_STATES", HOR_STATES);

	return fgm;
}

PyObject* python_create_alm(const struct mpc_ctl* data) {
	PyObject* alm = PyDict_New();

	PyObject* fgm = python_create_fgm(data);
	PyDict_SetItemString(alm, "fgm", fgm);

	python_add_matrix(
		alm, "l_0", data->alm->l_0,
		data->alm->HOR_MXCONSTRS, 1
	);

	python_add_matrix(
		alm, "zx_lb", data->alm->zx_lb,
		data->alm->HOR_MXCONSTRS, 1
	);

	python_add_matrix(
		alm, "zx_ub", data->alm->zx_ub,
		data->alm->HOR_MXCONSTRS, 1
	);

	PyObject* i_ex = PyLong_FromLong(*(data->alm->i_ex));
	PyDict_SetItemString(alm, "i_ex", i_ex);

	PyObject* mu = PyFloat_FromDouble(*(data->alm->mu));
	PyDict_SetItemString(alm, "mu", mu);

	python_add_matrix(
		alm, "E", data->alm->E,
		data->alm->HOR_MXCONSTRS, data->alm->HOR_INPUTS
	);

	python_add_matrix(
		alm, "Kx_Ai", data->alm->Kx_Ai,
		data->alm->HOR_MXCONSTRS, data->alm->STATES
	);

	python_add_matrix(
		alm, "e_lb", data->alm->e_lb,
		data->alm->HOR_MXCONSTRS, 1
	);

	python_add_matrix(
		alm, "e_ub", data->alm->e_ub,
		data->alm->HOR_MXCONSTRS, 1
	);

	PyObject* Linv = PyFloat_FromDouble(*(data->alm->Linv));
	PyDict_SetItemString(alm, "Linv", Linv);

	PyObject* STATES = PyFloat_FromDouble(data->alm->STATES);
	PyDict_SetItemString(alm, "STATES", STATES);

	PyObject* MXCONSTRS = PyFloat_FromDouble(data->alm->MXCONSTRS);
	PyDict_SetItemString(alm, "MXCONSTRS", MXCONSTRS);

	PyObject* HOR_INPUTS = PyFloat_FromDouble(data->alm->HOR_INPUTS);
	PyDict_SetItemString(alm, "HOR_INPUTS", HOR_INPUTS);

	PyObject* HOR_MXCONSTRS = PyFloat_FromDouble(data->alm->HOR_MXCONSTRS);
	PyDict_SetItemString(alm, "HOR_MXCONSTRS", HOR_MXCONSTRS);

	return alm;
}

PyObject* python_create_ctl(const struct mpc_ctl* data) {
	PyObject* ctl = PyDict_New();

	PyObject* conf = python_create_conf(data);
	PyDict_SetItemString(ctl, "conf", conf);

	PyObject* sys = python_create_sys(data);
	PyDict_SetItemString(ctl, "sys", sys);

	PyObject* wmx = python_create_wmx(data);
	PyDict_SetItemString(ctl, "wmx", wmx);

	PyObject* lqr = python_create_lqr(data);
	PyDict_SetItemString(ctl, "lqr", lqr);

	PyObject* qpx = python_create_qpx(data);
	PyDict_SetItemString(ctl, "qpx", qpx);

	PyObject* alm = python_create_alm(data);
	PyDict_SetItemString(ctl, "alm", alm);

	python_add_matrix(ctl, "u_opt", data->u_opt,
		data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(ctl, "l_opt", data->l_opt,
		data->alm->HOR_MXCONSTRS, 1
	);

	python_add_matrix(ctl, "x_trj", data->x_trj,
		data->alm->fgm->HOR_STATES+data->alm->STATES, 1
	);

	python_add_matrix(ctl, "u_ini", data->u_ini,
		data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(ctl, "l_ini", data->l_ini,
		data->alm->HOR_MXCONSTRS, 1
	);

	python_add_matrix(ctl, "u_ref", data->u_ref,
			data->alm->fgm->HOR_INPUTS, 1
	);

	python_add_matrix(ctl, "x_ref", data->x_ref,
		data->alm->fgm->HOR_STATES, 1
	);

	return ctl;
}
