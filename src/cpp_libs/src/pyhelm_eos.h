#ifndef PYHELMEOS_H
#define PYHELMEOS_H

#include "helm_eos.h"

typedef struct {
    PyObject_HEAD;
    t_helm_eos_table *helm_eos_table;
} pyHelmEos;

void pyHelmEosDealloc(PyObject* self);
PyObject* pyHelmEos_str(pyHelmEos* self);
PyObject* pyHelmEos_egiven(pyHelmEos* self, PyObject* args);
PyObject* pyHelmEos_pgiven(pyHelmEos* self, PyObject* args);
PyObject* pyHelmEos_tgiven(pyHelmEos* self, PyObject* args);
PyObject* pyHelmEos_tgivenfull(pyHelmEos* self, PyObject* args);
PyObject* pyHelmEos_egivenfull(pyHelmEos* self, PyObject* args);
PyObject* pyHelmEos_ptgivenfull(pyHelmEos* self, PyObject* args);
PyObject* PyGetHelmEosObject(t_helm_eos_table* helm_eos_table);
int pyConvertHelmEos(PyObject* object, t_helm_eos_table** helm_eos_table);
PyObject* loadhelm_eos(PyObject *self, PyObject *args);

// Python Module declarations
static PyMethodDef pyHelmEosMethods[] = {
    {"egiven", (getattrofunc)pyHelmEos_egiven, METH_VARARGS,
                           "egiven(rho, xnuc, e): returns temperature and pressure for given values of rho, xnuc and e"},
    {"pgiven", (getattrofunc)pyHelmEos_pgiven, METH_VARARGS,
                           "pgiven(rho, xnuc, p): returns temperature and energy for given values of rho, xnuc and p"},
    {"tgiven", (getattrofunc)pyHelmEos_tgiven, METH_VARARGS,
                           "tgiven(rho, xnuc, t): returns energy, dedT, pressure and sound speed for given values of rho, xnuc and t"},
    {"tgivenfull", (getattrofunc)pyHelmEos_tgivenfull, METH_VARARGS,
                           "tgivenfull(rho, xnuc, t): returns dictionary with all eos values for given values of rho, xnuc and t"},
    {"ptgivenfull", (getattrofunc)pyHelmEos_ptgivenfull, METH_VARARGS,
                           "ptgivenfull(p, xnuc, t [,rho]): returns dictionary with all eos values for given values of p, xnuc and t (initial guess for rho can be supplied)"},
    {"egivenfull", (getattrofunc)pyHelmEos_egivenfull, METH_VARARGS,""},
    {nullptr, nullptr, 0, nullptr}
};

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

static struct PyModuleDef moduledef_pyhelm_eos = {
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


#endif
