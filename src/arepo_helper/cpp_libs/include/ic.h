#ifndef AREPO_HELPER_LIBS_IC_H
#define AREPO_HELPER_LIBS_IC_H

#include <Python.h>
#include <arrayobject.h>
#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "helm_eos.h"
#include "utils.h"

struct paramsWD {
    double temp;
    double *xnuc;
    t_helm_eos_table *eos;
    double rho;
};

double rho_c_from_mtot_implementation(double mtot, double temp_c, t_helm_eos_table *eos, PyObject *xnuc_py);

double mtot_from_rho_c_implementation(double rho_c, double temp_c, t_helm_eos_table *eos, PyObject *xnuc_py, double offset = 0.0);

int create_wd_integrator(double r, const double *y, double *ydot, void *params);

PyObject *create_wd_implementation(t_helm_eos_table *eos,
                                   double rho_c,
                                   double temp_c,
                                   PyObject *xnuc_py,
                                   double tolerance);

PyObject *create_polytrope_implementation(t_helm_eos_table *eos,
                                          double n,
                                          double rho_c,
                                          PyObject *xnuc_py,
                                          double pres_c,
                                          double temp_c,
                                          double dr);

PyObject *create_wd_wdec_implementation(const char *wdec_dir, int nspecies, double gamma);

PyObject* create_wd(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject* create_polytrope(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject *create_wd_wdec(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject *rho_c_from_mtot(PyObject *self, PyObject *args);

PyObject *mtot_from_rho_c(PyObject *self, PyObject *args);

PyObject* _createWhiteDwarfHe(PyObject *self, PyObject *args);

// Python Module Definition
static PyMethodDef ic_methods[] = {
        {"create_wd", (PyCFunction) create_wd,METH_VARARGS | METH_KEYWORDS, ""},
        {"create_polytrope", (PyCFunction) create_polytrope, METH_VARARGS | METH_KEYWORDS, ""},
        {"create_wd_wdec", (PyCFunction) create_wd_wdec,METH_VARARGS | METH_KEYWORDS, ""},
        {"mtot_from_rho_c", (PyCFunction) mtot_from_rho_c, METH_VARARGS, ""},
        {"rho_c_from_mtot", (PyCFunction) rho_c_from_mtot, METH_VARARGS, ""},
        {"create_wdHe", _createWhiteDwarfHe, METH_VARARGS, "" },
        { nullptr, nullptr, 0, nullptr}
};

static struct PyModuleDef moduledef_ic = {
        PyModuleDef_HEAD_INIT,
        "ic",
        nullptr,
        -1,
        ic_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

#endif //AREPO_HELPER_LIBS_IC_H