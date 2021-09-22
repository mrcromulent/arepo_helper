#include <Python.h>
#include <arrayobject.h>

#ifndef AREPO_HELPER_LIBS_PYOPAL_EOS_H
#define AREPO_HELPER_LIBS_PYOPAL_EOS_H

#include "opal_eos.h"

int pyConvertOpalEos( PyObject* object, struct opal_eos_table** opal_eos_table );

#endif //AREPO_HELPER_LIBS_PYOPAL_EOS_H
