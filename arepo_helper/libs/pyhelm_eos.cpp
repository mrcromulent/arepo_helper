#include <Python.h>
#include <arrayobject.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pyhelm_eos.h"
#include "helm_eos.h"
#include "utils.h"


static PyMethodDef pyHelmEosMethods[] =
        {
                { "egiven", (getattrofunc)pyHelmEos_egiven, METH_VARARGS,
                        "egiven(rho, xnuc, e): returns temperature and pressure for given values of rho, xnuc and e" },
                { "pgiven", (getattrofunc)pyHelmEos_pgiven, METH_VARARGS,
                        "pgiven(rho, xnuc, p): returns temperature and energy for given values of rho, xnuc and p" },
                { "tgiven", (getattrofunc)pyHelmEos_tgiven, METH_VARARGS,
                        "tgiven(rho, xnuc, t): returns energy, dedT, pressure and sound speed for given values of rho, xnuc and t" },
                { "tgivenfull", (getattrofunc)pyHelmEos_tgivenfull, METH_VARARGS,
                        "tgivenfull(rho, xnuc, t): returns dictionary with all eos values for given values of rho, xnuc and t" },
                { "ptgivenfull", (getattrofunc)pyHelmEos_ptgivenfull, METH_VARARGS,
                        "ptgivenfull(p, xnuc, t [,rho]): returns dictionary with all eos values for given values of p, xnuc and t (initial guess for rho can be supplied)" },
                { nullptr, nullptr, 0, nullptr }
        };


static PyTypeObject pyHelmEosType = {
        PyVarObject_HEAD_INIT(nullptr, 0)
        "pyhelmeos",
        sizeof(pyHelmEos),
        0,
        pyHelmEosDealloc,
        0,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        (reprfunc)pyHelmEos_str,
        nullptr,
        nullptr,
        nullptr,
        Py_TPFLAGS_DEFAULT,
        "Helmholz EOS",
        nullptr,
        nullptr,
        nullptr,
        0,
        nullptr,
        nullptr,
        pyHelmEosMethods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        0,
        nullptr,
        nullptr,
        nullptr,
};

void pyHelmEosDealloc( PyObject* self ) {
    eos_deinit( ((pyHelmEos*)self)->helm_eos_table );
    PyObject_Del( self );
}

PyObject* pyHelmEos_str( pyHelmEos* self ) {
    return PyUnicode_FromString( "helm_eos" );
}

PyObject* pyHelmEos_egiven( pyHelmEos* self, PyObject* args ) {
    PyArrayObject *pyXnuc;
    double rho, e, temp;
    struct eos_result res{};

    if ( !PyArg_ParseTuple( args, "dO!d:helmeos.egiven(rho, xnuc, e)", &rho, &PyArray_Type, &pyXnuc, &e ) ) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        printf( "%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies );
        PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
        return nullptr;
    }

    temp = -1.;
    eos_calc_egiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), e, &temp, &res );

    return Py_BuildValue( "dd", temp, res.p.v );
}

PyObject* pyHelmEos_pgiven( pyHelmEos* self, PyObject* args ) {
    PyArrayObject *pyXnuc;
    double rho, p, temp;
    struct eos_result res{};

    if ( !PyArg_ParseTuple( args, "dO!d:helmeos.pgiven(rho, xnuc, p)", &rho, &PyArray_Type, &pyXnuc, &p ) ) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        printf( "%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies );
        PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
        return nullptr;
    }

    temp = -1.;
    eos_calc_pgiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), p, &temp, &res );

    return Py_BuildValue( "dd", temp, res.e.v );
}

PyObject* pyHelmEos_tgiven( pyHelmEos* self, PyObject* args ) {
    PyArrayObject *pyXnuc;
    double rho, temp;
    struct eos_result res{};

    if ( !PyArg_ParseTuple( args, "dO!d:helmeos.tgiven(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp ) ) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
        return nullptr;
    }

    eos_calc_tgiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &res );

    return Py_BuildValue( "dddd", res.e.v, res.e.dtemp, res.p.v, res.sound );
}

PyObject* pyHelmEos_tgivenfull( pyHelmEos* self, PyObject* args ) {
    PyArrayObject *pyXnuc;
    double rho, temp;
    struct eos_result res{};

    if ( !PyArg_ParseTuple( args, "dO!d:helmeos.tgivenfull(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp ) ) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
        return nullptr;
    }

    eos_calc_tgiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &res );

    PyObject* dict = PyDict_New();
    PyDict_SetStolenItem( dict, "e", (PyObject*)PyFloat_FromDouble( res.e.v ) );
    PyDict_SetStolenItem( dict, "dedT", (PyObject*)PyFloat_FromDouble( res.e.dtemp ) );
    PyDict_SetStolenItem( dict, "T", (PyObject*)PyFloat_FromDouble( res.temp ) );
    PyDict_SetStolenItem( dict, "p", (PyObject*)PyFloat_FromDouble( res.p.v ) );
    PyDict_SetStolenItem( dict, "dpdT", (PyObject*)PyFloat_FromDouble( res.p.dtemp ) );
    PyDict_SetStolenItem( dict, "csnd", (PyObject*)PyFloat_FromDouble( res.sound ) );
    PyDict_SetStolenItem( dict, "gamma1", (PyObject*)PyFloat_FromDouble( res.gamma_1 ) );
    PyDict_SetStolenItem( dict, "gamma2", (PyObject*)PyFloat_FromDouble( res.gamma_2 ) );
    PyDict_SetStolenItem( dict, "gamma3", (PyObject*)PyFloat_FromDouble( res.gamma_3 ) );
    PyDict_SetStolenItem( dict, "cv", (PyObject*)PyFloat_FromDouble( res.cv ) );
    PyDict_SetStolenItem( dict, "cp", (PyObject*)PyFloat_FromDouble( res.cp ) );

    return dict;
}

PyObject* pyHelmEos_ptgivenfull( pyHelmEos* self, PyObject* args ) {
    PyArrayObject *pyXnuc;
    double p, temp, rho;
    struct eos_result res{};

    rho = -1;
    if ( !PyArg_ParseTuple( args, "dO!d|d:helmeos.ptgivenfull(p, xnuc, temp, [rho])", &p, &PyArray_Type, &pyXnuc, &temp, &rho) ) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
        PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
        return nullptr;
    }

    eos_calc_ptgiven(self->helm_eos_table, p, (double*)PyArray_DATA(pyXnuc), temp, &rho, &res);

    PyObject* dict = PyDict_New();
    PyDict_SetStolenItem( dict, "rho", (PyObject*)PyFloat_FromDouble(rho));
    PyDict_SetStolenItem( dict, "e", (PyObject*)PyFloat_FromDouble( res.e.v ) );
    PyDict_SetStolenItem( dict, "dedT", (PyObject*)PyFloat_FromDouble( res.e.dtemp ) );
    PyDict_SetStolenItem( dict, "T", (PyObject*)PyFloat_FromDouble( res.temp ) );
    PyDict_SetStolenItem( dict, "p", (PyObject*)PyFloat_FromDouble( res.p.v ) );
    PyDict_SetStolenItem( dict, "dpdT", (PyObject*)PyFloat_FromDouble( res.p.dtemp ) );
    PyDict_SetStolenItem( dict, "csnd", (PyObject*)PyFloat_FromDouble( res.sound ) );
    PyDict_SetStolenItem( dict, "gamma1", (PyObject*)PyFloat_FromDouble( res.gamma_1 ) );
    PyDict_SetStolenItem( dict, "gamma2", (PyObject*)PyFloat_FromDouble( res.gamma_2 ) );
    PyDict_SetStolenItem( dict, "gamma3", (PyObject*)PyFloat_FromDouble( res.gamma_3 ) );
    PyDict_SetStolenItem( dict, "cv", (PyObject*)PyFloat_FromDouble( res.cv ) );
    PyDict_SetStolenItem( dict, "cp", (PyObject*)PyFloat_FromDouble( res.cp ) );

    return dict;
}

PyObject* PyGetHelmEosObject( t_helm_eos_table* helm_eos_table ) {
    pyHelmEos* eos;

    if (helm_eos_table) {
        eos = (pyHelmEos*)PyObject_New( pyHelmEos, &pyHelmEosType );
        eos->helm_eos_table = helm_eos_table;
        return (PyObject*)eos;
    } else {
        Py_RETURN_NONE;
    }
}

int pyConvertHelmEos( PyObject* object, t_helm_eos_table** helm_eos_table ) {
    if (strcmp(object->ob_type->tp_name, "pyhelmeos" ) != 0) {
        PyErr_BadArgument();
        return 0;
    }

    *helm_eos_table = ((pyHelmEos*)object)->helm_eos_table;
    return 1;
}

/* pyEos module */

PyObject* loadhelm_eos(PyObject *self, PyObject *args) {
    char *datafile, *speciesfile;
    t_helm_eos_table* helm_eos_table;

    int useCoulombCorrections = 0;
    int verbose = 0;
    if (!PyArg_ParseTuple( args, "ss|ii:loadhelm_eos( datafile, speciesfile, [useCoulombCorrections,verbose] )", &datafile, &speciesfile, &useCoulombCorrections, &verbose )) {
        return nullptr;
    }

    helm_eos_table = (t_helm_eos_table*)malloc( sizeof( t_helm_eos_table ) );
    eos_init( helm_eos_table, datafile, speciesfile, useCoulombCorrections, verbose );

    return PyGetHelmEosObject( helm_eos_table );
}

// Python Module declarations
static PyMethodDef pyhelmeosMethods[] = {
        {"loadhelm_eos",
                loadhelm_eos,
        METH_VARARGS,
        "loadhelm_eos(datafile, speciesfile, [useCoulombCorrections, verbose]).\n Load Helmholz EOS from specified datafile and speciesfile."
        },
        { nullptr,
          nullptr,
          0,
          nullptr
        }
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pyhelm_eos",
        nullptr,
        -1,
        pyhelmeosMethods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_pyhelm_eos(void) {

    import_array();
    pyHelmEosType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyHelmEosType) < 0) {
        return nullptr;
    }

    return PyModule_Create(&moduledef);
}
