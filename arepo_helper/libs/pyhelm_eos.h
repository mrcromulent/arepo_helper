#ifndef PYHELMEOS_H
#define PYHELMEOS_H

#include "helm_eos.h"

int pyConvertHelmEos( PyObject* object, t_helm_eos_table** helm_eos_table );

#endif
