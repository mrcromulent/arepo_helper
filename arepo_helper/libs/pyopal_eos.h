#ifndef PYOPALEOS_H
#define PYOPALEOS_H

#include "opal_eos.h"

int pyConvertOpalEos( PyObject* object, struct opal_eos_table** opal_eos_table );

#endif
