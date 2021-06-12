#ifndef AREPO_HELPER_LIBS_IC_H
#define AREPO_HELPER_LIBS_IC_H

#include <Python.h>
#include <arrayobject.h>
#include "helm_eos.h"
#include <algorithm>


struct paramsWD {
    double temp;
    double *xnuc;
    t_helm_eos_table *eos;
    double rho;
};

int createWDIntegrator( double r, const double y[], double f[], void *params);
PyObject* _createWhiteDwarf(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject* _createPolytrope(PyObject *self, PyObject *args, PyObject *kwargs);


#endif //AREPO_HELPER_LIBS_IC_H