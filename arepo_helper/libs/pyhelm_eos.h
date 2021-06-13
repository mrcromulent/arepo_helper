#ifndef PYHELMEOS_H
#define PYHELMEOS_H

#include "helm_eos.h"

typedef struct {
    PyObject_HEAD;
    t_helm_eos_table *helm_eos_table;
} pyHelmEos;

/* forward declarations */
void pyHelmEosDealloc( PyObject* self );
static PyObject* pyHelmEos_str( pyHelmEos* self );
PyObject* pyHelmEos_egiven( pyHelmEos* self, PyObject* args );
PyObject* pyHelmEos_pgiven( pyHelmEos* self, PyObject* args );
PyObject* pyHelmEos_tgiven( pyHelmEos* self, PyObject* args );
PyObject* pyHelmEos_tgivenfull( pyHelmEos* self, PyObject* args );
PyObject* pyHelmEos_ptgivenfull( pyHelmEos* self, PyObject* args );
int pyConvertHelmEos( PyObject* object, t_helm_eos_table** helm_eos_table );

#endif
