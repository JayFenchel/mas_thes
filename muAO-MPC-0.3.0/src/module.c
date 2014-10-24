#include "cmpc.h"

#define PY_ARRAY_UNIQUE_SYMBOL _muaompc_cmpc_numpy
#include <numpy/npy_common.h>
#include <numpy/arrayobject.h>

struct mpc_vla {
	real_t* HoL;
	real_t* GoL;
	real_t* Bh_T;
	real_t* E;
	real_t* Kx_Ai;
	real_t* u_lb;
	real_t* u_ub;
	real_t* e_lb;
	real_t* e_ub;
	real_t nu;
	real_t mu;
	real_t Linv;
	real_t* gxoL;
	real_t* zx_lb;
	real_t* zx_ub;
	real_t* groL;
	real_t* Q;
	real_t* R;
	real_t* P;
	real_t* K;
	real_t* Ad;
	real_t* Bd;
	real_t dt;
	real_t* u_opt;
	real_t* l_opt;
	real_t* x_trj;
	real_t* u_ini;
	real_t* l_ini;
	real_t* u_ref;
	real_t* x_ref;

	uint32_t HOR_INPUTS;
	uint32_t HOR_MXCONSTRS;
	uint32_t HOR_STATES;
	uint32_t HOR;
	uint32_t INPUTS;
	uint32_t MXCONSTRS;
	uint32_t STATES;

	struct mpc_conf conf;
	struct mpc_qpx qpx;
	struct mpc_sys sys;
	struct mpc_wmx wmx;
	struct mpc_lqr lqr;
	struct mpc_fgm fgm;
	struct mpc_alm alm;
	struct mpc_ctl ctl;

	real_t* x;
} vla;

real_t* vla_alloc(uint32_t size) {
	return (real_t*)calloc(size, sizeof(real_t));
}

void vla_free(real_t* p) {
	free((void*)p);
}

void vla_init(struct mpc_vla* vla, PyObject* ctl) {
	PyObject* alm = PyDict_GetItemString(ctl, "alm");
	PyObject* fgm = PyDict_GetItemString(alm, "fgm");
	PyObject* sys = PyDict_GetItemString(ctl, "sys");
	PyObject* wmx = PyDict_GetItemString(ctl, "wmx");
	PyObject* lqr = PyDict_GetItemString(ctl, "lqr");
	PyObject* conf = PyDict_GetItemString(ctl, "conf");
	PyObject* qpx = PyDict_GetItemString(ctl, "qpx");

	PyObject* in_iter = PyDict_GetItemString(conf, "in_iter");
	vla->conf.in_iter = PyLong_AsLong(in_iter);

	PyObject* ex_iter = PyDict_GetItemString(conf, "ex_iter");
	vla->conf.ex_iter = PyLong_AsLong(ex_iter);

	PyObject* warmstart = PyDict_GetItemString(conf, "warmstart");
	vla->conf.warmstart = PyLong_AsLong(warmstart);

	PyObject* HOR_INPUTS = PyDict_GetItemString(alm, "HOR_INPUTS");
	vla->HOR_INPUTS = PyLong_AsLong(HOR_INPUTS);

	PyObject* HOR_MXCONSTRS = PyDict_GetItemString(alm, "HOR_MXCONSTRS");
	vla->HOR_MXCONSTRS = PyLong_AsLong(HOR_MXCONSTRS);

	PyObject* HOR_STATES = PyDict_GetItemString(fgm, "HOR_STATES");
	vla->HOR_STATES = PyLong_AsLong(HOR_STATES);

	PyObject* HOR = PyDict_GetItemString(fgm, "HOR");
	vla->HOR = PyLong_AsLong(HOR);

	PyObject* INPUTS = PyDict_GetItemString(fgm, "INPUTS");
	vla->INPUTS = PyLong_AsLong(INPUTS);

	PyObject* MXCONSTRS = PyDict_GetItemString(alm, "MXCONSTRS");
	vla->MXCONSTRS = PyLong_AsLong(MXCONSTRS);

	PyObject* STATES = PyDict_GetItemString(alm, "STATES");
	vla->STATES = PyLong_AsLong(STATES);

	vla->HoL = vla_alloc(vla->HOR_INPUTS * vla->HOR_INPUTS);
	python_read_matrix(
		qpx, "HoL", vla->HoL,
		vla->HOR_INPUTS, vla->HOR_INPUTS
	);

	vla->GoL = vla_alloc(vla->HOR_INPUTS * vla->STATES);
	python_read_matrix(
		fgm, "GoL", vla->GoL,
		vla->HOR_INPUTS, vla->STATES
	);

	vla->Bh_T = vla_alloc(vla->HOR_INPUTS * vla->HOR_STATES);
	python_read_matrix(
		fgm, "Bh_T", vla->Bh_T,
		vla->HOR_INPUTS, vla->HOR_STATES
	);

	vla->E = vla_alloc(vla->HOR_MXCONSTRS * vla->HOR_INPUTS);
	python_read_matrix(
		alm, "E", vla->E,
		vla->HOR_MXCONSTRS, vla->HOR_INPUTS
	);

	vla->Kx_Ai = vla_alloc(vla->HOR_MXCONSTRS * vla->STATES);
	python_read_matrix(
		alm, "Kx_Ai", vla->Kx_Ai,
		vla->HOR_MXCONSTRS, vla->STATES
	);

	vla->u_lb = vla_alloc(vla->HOR_INPUTS);
	python_read_matrix(
		fgm, "u_lb", vla->u_lb,
		vla->HOR_INPUTS, 1
	);

	vla->u_ub = vla_alloc(vla->HOR_INPUTS);
	python_read_matrix(
		fgm, "u_ub", vla->u_ub,
		vla->HOR_INPUTS, 1
	);

	vla->e_lb = vla_alloc(vla->HOR_MXCONSTRS);
	python_read_matrix(
		alm, "e_lb", vla->e_lb,
		vla->HOR_MXCONSTRS, 1
	);

	vla->e_ub = vla_alloc(vla->HOR_MXCONSTRS);
	python_read_matrix(
		alm, "e_ub", vla->e_ub,
		vla->HOR_MXCONSTRS, 1
	);

	PyObject* nu = PyDict_GetItemString(fgm, "nu");
	vla->nu = PyFloat_AsDouble(nu);

	PyObject* mu = PyDict_GetItemString(alm, "mu");
	vla->mu = PyFloat_AsDouble(mu);

	PyObject* Linv = PyDict_GetItemString(alm, "Linv");
	vla->Linv = PyFloat_AsDouble(Linv);

	vla->gxoL = vla_alloc(vla->HOR_INPUTS);
	vla->zx_lb = vla_alloc(vla->HOR_MXCONSTRS);
	vla->zx_ub = vla_alloc(vla->HOR_MXCONSTRS);
	vla->groL = vla_alloc(vla->HOR_INPUTS);

	vla->Q = vla_alloc(vla->STATES * vla->STATES);
	python_read_matrix(
		wmx, "Q", vla->Q,
		vla->STATES, vla->STATES
	);

	vla->R = vla_alloc(vla->INPUTS * vla->INPUTS);
	python_read_matrix(
		wmx, "R", vla->R,
		vla->INPUTS, vla->INPUTS
	);

	vla->P = vla_alloc(vla->STATES * vla->STATES);
	python_read_matrix(
		wmx, "P", vla->P,
		vla->STATES, vla->STATES
	);

	vla->K = vla_alloc(vla->INPUTS * vla->STATES);
	python_read_matrix(
		lqr, "K", vla->K,
		vla->INPUTS, vla->STATES
	);

	vla->P = vla_alloc(vla->STATES * vla->STATES);
	python_read_matrix(
		lqr, "P", vla->P,
		vla->STATES, vla->STATES
	);

	vla->Ad = vla_alloc(vla->STATES * vla->STATES);
	python_read_matrix(
		sys, "Ad", vla->Ad,
		vla->STATES, vla->STATES
	);

	vla->Bd = vla_alloc(vla->STATES * vla->INPUTS);
	python_read_matrix(
		sys, "Bd", vla->Bd,
		vla->STATES, vla->INPUTS
	);

	PyObject* dt = PyDict_GetItemString(sys, "dt");
	vla->dt = PyFloat_AsDouble(dt);

	vla->u_opt = vla_alloc(vla->HOR_INPUTS);
	vla->l_opt = vla_alloc(vla->HOR_MXCONSTRS);
	vla->x_trj = vla_alloc(vla->HOR_STATES+vla->STATES);

	vla->u_ini = vla_alloc(vla->HOR_INPUTS);
	python_read_matrix(
		ctl, "u_ini", vla->u_ini,
		vla->HOR_INPUTS, 1
	);

	vla->l_ini = vla_alloc(vla->HOR_MXCONSTRS);
	python_read_matrix(
		ctl, "l_ini", vla->l_ini,
		vla->HOR_MXCONSTRS, 1
	);

	PyObject* u_ref = PyDict_GetItemString(ctl, "u_ref");
	vla->u_ref = vla_alloc(vla->HOR_INPUTS);
	python_read_matrix(
		ctl, "u_ref", vla->u_ref,
			vla->HOR_INPUTS, 1
	);

	PyObject* x_ref = PyDict_GetItemString(ctl, "x_ref");
	vla->x_ref = vla_alloc(vla->HOR_STATES);
	python_read_matrix(
		ctl, "x_ref", vla->x_ref,
		vla->HOR_STATES, 1
	);

	vla->x = vla_alloc(vla->STATES);

	vla->qpx.HoL = vla->HoL;
	vla->qpx.gxoL = vla->gxoL;
	vla->qpx.E = vla->E;
	vla->qpx.u_lb = vla->u_lb;
	vla->qpx.u_ub = vla->u_ub;
	vla->qpx.zx_lb = vla->zx_lb;
	vla->qpx.zx_ub = vla->zx_ub;
	vla->qpx.HOR_INPUTS = vla->HOR_INPUTS;
	vla->qpx.HOR_MXCONSTRS = vla->HOR_MXCONSTRS;

	vla->sys.Ad = vla->Ad;
	vla->sys.Bd = vla->Bd;
	vla->sys.dt = &(vla->dt);

	vla->wmx.Q = vla->Q;
	vla->wmx.R = vla->R;
	vla->wmx.P = vla->P;

	vla->lqr.K = vla->K;
	vla->lqr.P = vla->P;

	vla->fgm.u_0 = vla->u_ini;
	vla->fgm.gxoL = vla->gxoL;
	vla->fgm.groL = vla->groL;
	vla->fgm.j_in = &(vla->conf.in_iter);
	vla->fgm.HoL = vla->HoL;
	vla->fgm.GoL = vla->GoL;
	vla->fgm.Bh_T = vla->Bh_T;
	vla->fgm.u_lb = vla->u_lb;
	vla->fgm.u_ub = vla->u_ub;
	vla->fgm.nu = &(vla->nu);
	vla->fgm.HOR = vla->HOR;
	vla->fgm.STATES = vla->STATES;
	vla->fgm.INPUTS = vla->INPUTS;
	vla->fgm.HOR_INPUTS = vla->HOR_INPUTS;
	vla->fgm.HOR_STATES = vla->HOR_STATES;

	vla->alm.fgm = &(vla->fgm);
	vla->alm.l_0 = vla->l_ini;
	vla->alm.zx_lb = vla->zx_lb;
	vla->alm.zx_ub = vla->zx_ub;
	vla->alm.i_ex = &(vla->conf.ex_iter);
	vla->alm.mu = &(vla->mu);
	vla->alm.E = vla->E;
	vla->alm.Kx_Ai = vla->Kx_Ai;
	vla->alm.e_lb = vla->e_lb;
	vla->alm.e_ub = vla->e_ub;
	vla->alm.Linv = &(vla->Linv);
	vla->alm.STATES = vla->STATES;
	vla->alm.MXCONSTRS = vla->MXCONSTRS;
	vla->alm.HOR_INPUTS = vla->HOR_INPUTS;
	vla->alm.HOR_MXCONSTRS = vla->HOR_MXCONSTRS;

	vla->ctl.conf = &(vla->conf);
	vla->ctl.qpx = &(vla->qpx);
	vla->ctl.sys = &(vla->sys);
	vla->ctl.wmx = &(vla->wmx);
	vla->ctl.lqr = &(vla->lqr);
	vla->ctl.alm = &(vla->alm);
	vla->ctl.u_opt = vla->u_opt;
	vla->ctl.l_opt = vla->l_opt;
	vla->ctl.x_trj = vla->x_trj;
	vla->ctl.u_ini = vla->u_ini;
	vla->ctl.l_ini = vla->l_ini;
	vla->ctl.u_ref = vla->u_ref;
	vla->ctl.x_ref = vla->x_ref;
}

void vla_cleanup(struct mpc_vla* vla) {
	vla_free(vla->HoL);
	vla_free(vla->GoL);
	vla_free(vla->Bh_T);
	vla_free(vla->E);
	vla_free(vla->Kx_Ai);
	vla_free(vla->u_lb);
	vla_free(vla->u_ub);
	vla_free(vla->e_lb);
	vla_free(vla->e_ub);
	vla_free(vla->gxoL);
	vla_free(vla->zx_lb);
	vla_free(vla->zx_ub);
	vla_free(vla->groL);
	vla_free(vla->Q);
	vla_free(vla->R);
	vla_free(vla->P);
	vla_free(vla->Ad);
	vla_free(vla->Bd);
	vla_free(vla->u_opt);
	vla_free(vla->l_opt);
	vla_free(vla->x_trj);
	vla_free(vla->u_ini);
	vla_free(vla->l_ini);
	vla_free(vla->u_ref);
	vla_free(vla->x_ref);
	vla_free(vla->x);
}

static PyObject* vla_form_qp(PyObject* self, PyObject* args) {
	PyObject* structure;
	PyArrayObject* x;

	if (!PyArg_ParseTuple(args, "O!O!", &PyDict_Type, &structure, &PyArray_Type, &x)) {
		Py_RETURN_NONE;
	}

	if (0 == x) {
		Py_RETURN_NONE;
	}

	if (0 == structure) {
		Py_RETURN_NONE;
	}

	vla_init(&vla, structure);
	memcpy(vla.x, PyArray_DATA(x), sizeof(float64_t) * vla.STATES);
	mpc_ctl_form_qp(&(vla.ctl), vla.x);
	PyObject* ret = python_create_ctl(&(vla.ctl));
	vla_cleanup(&vla);
	return ret;
}

static PyObject* vla_solve_problem(PyObject* self, PyObject* args) {
	PyObject* structure;
	PyArrayObject* x;

	if (!PyArg_ParseTuple(args, "O!O!", &PyDict_Type, &structure, &PyArray_Type, &x)) {
		Py_RETURN_NONE;
	}

	if (0 == x) {
		Py_RETURN_NONE;
	}

	if (0 == structure) {
		Py_RETURN_NONE;
	}

	vla_init(&vla, structure);
	memcpy(vla.x, PyArray_DATA(x), sizeof(float64_t) * vla.STATES);
	mpc_ctl_solve_problem(&(vla.ctl), vla.x);
	PyObject* ret = python_create_ctl(&(vla.ctl));
	vla_cleanup(&vla);
	return ret;
}

static PyObject* vla_get_data(PyObject* self, PyObject* args) {
	PyObject* structure;

	if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &structure)) {
		Py_RETURN_NONE;
	}

	if (0 == structure) {
		Py_RETURN_NONE;
	}

	vla_init(&vla, structure);
	PyObject* ret = python_create_ctl(&(vla.ctl));
	vla_cleanup(&vla);
	return ret;
}

static PyMethodDef methods[] = {
	{"vla_get_data", vla_get_data, METH_VARARGS, ""},
	{"vla_solve_problem", vla_solve_problem, METH_VARARGS, ""},
	{"vla_form_qp", vla_form_qp, METH_VARARGS, ""},
	{0, 0, 0, 0}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef module = {
	PyModuleDef_HEAD_INIT,
	"_cmpc",
	0,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__cmpc() {
	import_array();
	return PyModule_Create(&module);
}

#elif PY_MAJOR_VERSION == 2
/* the following allows to compile for Python2,
 * however, the ctl will not be available inside Python
 * more tests are required */

PyMODINIT_FUNC

init_cmpc(void)
{
	import_array();
     (void) Py_InitModule("_cmpc", methods);
}

#include "Python.h"
int main()
{
 Py_Initialize();
 Py_InitModule( "_cmpc", methods );
 Py_Finalize();
 return 0;
}
#endif
