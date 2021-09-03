#include <Python.h>
#include <arrayobject.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pyhelm_eos.h"
#include "helm_eos.h"
#include "utils.h"
#include "const.h"


void pyHelmEosDealloc(PyObject* self) {
    eos_deinit(((pyHelmEos*)self)->helm_eos_table);
    PyObject_Del(self);
}

PyObject* pyHelmEos_str(pyHelmEos* self) {
    return PyUnicode_FromString("helm_eos");
}

PyObject* pyHelmEos_egiven(pyHelmEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, e, temp;
    struct eos_result res{};

    if (!PyArg_ParseTuple(args, "dO!d:helmeos.egiven(rho, xnuc, e)", &rho, &PyArray_Type, &pyXnuc, &e)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        printf("%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies);
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    temp = -1.;
    eos_calc_egiven(self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), e, &temp, &res);

    return Py_BuildValue("dd", temp, res.p.v);
}

PyObject* pyHelmEos_pgiven(pyHelmEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, p, temp;
    struct eos_result res{};

    if (!PyArg_ParseTuple(args, "dO!d:helmeos.pgiven(rho, xnuc, p)", &rho, &PyArray_Type, &pyXnuc, &p)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        printf("%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies);
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    temp = -1.;
    eos_calc_pgiven(self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), p, &temp, &res);

    return Py_BuildValue("dd", temp, res.e.v);
}

PyObject* pyHelmEos_tgiven(pyHelmEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, temp;
    struct eos_result res{};

    if (!PyArg_ParseTuple(args, "dO!d:helmeos.tgiven(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    eos_calc_tgiven(self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &res);

    return Py_BuildValue("dddd", res.e.v, res.e.dtemp, res.p.v, res.sound);
}

PyObject* pyHelmEos_tgivenfull(pyHelmEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, temp;
    struct eos_result res{};

    if (!PyArg_ParseTuple(args, "dO!d:helmeos.tgivenfull(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    eos_calc_tgiven(self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &res);

    PyObject* dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject*)PyFloat_FromDouble(res.e.v));
    PyDict_SetStolenItem(dict, f[N::DEDT], (PyObject*)PyFloat_FromDouble(res.e.dtemp));
    PyDict_SetStolenItem(dict, f[N::TEMPERATURE], (PyObject*)PyFloat_FromDouble(res.temp));
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject*)PyFloat_FromDouble(res.p.v));
    PyDict_SetStolenItem(dict, f[N::DPDT], (PyObject*)PyFloat_FromDouble(res.p.dtemp));
    PyDict_SetStolenItem(dict, f[N::SOUNDSPEED], (PyObject*)PyFloat_FromDouble(res.sound));
    PyDict_SetStolenItem(dict, f[N::GAMMA1], (PyObject*)PyFloat_FromDouble(res.gamma_1));
    PyDict_SetStolenItem(dict, f[N::GAMMA2], (PyObject*)PyFloat_FromDouble(res.gamma_2));
    PyDict_SetStolenItem(dict, f[N::GAMMA3], (PyObject*)PyFloat_FromDouble(res.gamma_3));
    PyDict_SetStolenItem(dict, f[N::CV], (PyObject*)PyFloat_FromDouble(res.cv));
    PyDict_SetStolenItem(dict, f[N::CP], (PyObject*)PyFloat_FromDouble(res.cp));

    return dict;
}

PyObject* pyHelmEos_ptgivenfull(pyHelmEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double p, temp, rho;
    struct eos_result res{};

    rho = -1;
    if (!PyArg_ParseTuple(args, "dO!d|d:helmeos.ptgivenfull(p, xnuc, temp, [rho])", &p, &PyArray_Type, &pyXnuc, &temp, &rho)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    eos_calc_ptgiven(self->helm_eos_table, p, (double*)PyArray_DATA(pyXnuc), temp, &rho, &res);

    PyObject* dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::DENSITY], (PyObject*)PyFloat_FromDouble(rho));
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject*)PyFloat_FromDouble(res.e.v));
    PyDict_SetStolenItem(dict, f[N::DEDT], (PyObject*)PyFloat_FromDouble(res.e.dtemp));
    PyDict_SetStolenItem(dict, f[N::TEMPERATURE], (PyObject*)PyFloat_FromDouble(res.temp));
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject*)PyFloat_FromDouble(res.p.v));
    PyDict_SetStolenItem(dict, f[N::DPDT], (PyObject*)PyFloat_FromDouble(res.p.dtemp));
    PyDict_SetStolenItem(dict, f[N::SOUNDSPEED], (PyObject*)PyFloat_FromDouble(res.sound));
    PyDict_SetStolenItem(dict, f[N::GAMMA1], (PyObject*)PyFloat_FromDouble(res.gamma_1));
    PyDict_SetStolenItem(dict, f[N::GAMMA2], (PyObject*)PyFloat_FromDouble(res.gamma_2));
    PyDict_SetStolenItem(dict, f[N::GAMMA3], (PyObject*)PyFloat_FromDouble(res.gamma_3));
    PyDict_SetStolenItem(dict, f[N::CV], (PyObject*)PyFloat_FromDouble(res.cv));
    PyDict_SetStolenItem(dict, f[N::CP], (PyObject*)PyFloat_FromDouble(res.cp));

    return dict;
}

PyObject* PyGetHelmEosObject(t_helm_eos_table* helm_eos_table) {
    pyHelmEos* eos;

    if (helm_eos_table) {
        eos = (pyHelmEos*)PyObject_New(pyHelmEos, &pyHelmEosType);
        eos->helm_eos_table = helm_eos_table;
        return (PyObject*)eos;
    } else {
        Py_RETURN_NONE;
    }
}

int pyConvertHelmEos(PyObject* object, t_helm_eos_table** helm_eos_table) {
    if (strcmp(object->ob_type->tp_name, "pyhelmeos") != 0) {
        PyErr_BadArgument();
        return 0;
    }

    *helm_eos_table = ((pyHelmEos*)object)->helm_eos_table;
    return 1;
}

PyObject* loadhelm_eos(PyObject *self, PyObject *args) {
    char *datafile, *speciesfile;
    t_helm_eos_table* helm_eos_table;

    int use_coulomb_corrections = 0;
    int verbose = 0;
    if (!PyArg_ParseTuple(args, "ss|ii:loadhelm_eos(datafile, speciesfile, [use_coulomb_corrections, verbose])",
                          &datafile, &speciesfile, &use_coulomb_corrections, &verbose)) {
        return nullptr;
    }

    helm_eos_table = (t_helm_eos_table*)malloc(sizeof(t_helm_eos_table));
    eos_init(helm_eos_table, datafile, speciesfile, use_coulomb_corrections, verbose);

    return PyGetHelmEosObject(helm_eos_table);
}

PyMODINIT_FUNC PyInit_pyhelm_eos(void) {

    import_array()
    pyHelmEosType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyHelmEosType) < 0) {
        return nullptr;
    }

    return PyModule_Create(&moduledef_pyhelm_eos);
}
