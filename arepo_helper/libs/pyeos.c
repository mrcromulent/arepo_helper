#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pyeos.h"
#include "eos.h"

typedef struct {
  PyObject_HEAD;
  t_eos_table *eos_table;
} pyEos;

/* forward declarations */
void pyEosDealloc( PyObject* self );
static PyObject* pyEos_str( pyEos* self );
PyObject* pyEos_egiven( pyEos* self, PyObject* args );
PyObject* pyEos_tgiven( pyEos* self, PyObject* args );

static PyMethodDef pyEosMethods[] =
{
{ "egiven", (getattrofunc)pyEos_egiven, METH_VARARGS, 0 },
{ "tgiven", (getattrofunc)pyEos_tgiven, METH_VARARGS, 0 },
{ NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3
static PyTypeObject pyEosType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "pyeos",               /* tp_name */
    sizeof(pyEos),         /* tp_basicsize */
    0,                         /* tp_itemsize */
    pyEosDealloc,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    (reprfunc)pyEos_str,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "old LEAFS EOS",            /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    pyEosMethods,          /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    0,                         /* tp_new */
};
#else
PyObject* pyEos_getAttr( pyEos* self, char* name ) {
  return Py_FindMethod( pyEosMethods, (PyObject*)self, name );
}

static PyTypeObject pyEosType =
{
PyObject_HEAD_INIT( &PyType_Type )
0,
"pyeos",
sizeof( pyEosType ),
0,
pyEosDealloc,
0,
( getattrfunc ) pyEos_getAttr,
0,0,0,0,0,0,0,0,
( reprfunc ) pyEos_str,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
};
#endif

void pyEosDealloc( PyObject* self ) {
  eos_deinit( ((pyEos*)self)->eos_table );
  PyObject_Del( self );
}

PyObject* pyEos_str( pyEos* self ) {
  return PyUnicode_FromString( self->eos_table->datafile );
}

PyObject* pyEos_egiven( pyEos* self, PyObject* args ) {
  PyArrayObject *pyXnuc;
  double rho, e, temp, p, dpdr;

  if ( !PyArg_ParseTuple( args, "dO!d:eos.egiven(rho, xnuc, e)", &rho, &PyArray_Type, &pyXnuc, &e ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->eos_table->nspecies) {
    printf( "%d %d\n", PyArray_NDIM(pyXnuc), self->eos_table->nspecies );
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
  }
  
  eos_calc_egiven( self->eos_table, rho, (double*)PyArray_DATA(pyXnuc), e, &temp, &p, &dpdr );
  
  return Py_BuildValue( "ddd", temp, p, dpdr );
}

PyObject* pyEos_tgiven( pyEos* self, PyObject* args ) {
  PyArrayObject *pyXnuc;
  double rho, temp, e, dedt, p;

  if ( !PyArg_ParseTuple( args, "dO!d:eos.tgiven(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->eos_table->nspecies) {
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
  }
  
  eos_calc_tgiven( self->eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &e, &dedt, &p );
  
  return Py_BuildValue( "ddd", e, dedt, p );
}

PyObject* PyGetEosObject( t_eos_table* eos_table ) {
  pyEos* eos;
  
  if (eos_table) {
    eos = (pyEos*)PyObject_New( pyEos, &pyEosType );
    eos->eos_table = eos_table;
    return (PyObject*)eos;
  } else {
    Py_RETURN_NONE;
  }
}

int pyConvertEos( PyObject* object, t_eos_table** eos_table ) {
  if (strcmp(object->ob_type->tp_name, "pyeos" )) {
    PyErr_BadArgument();
    return 0;
  }
  
  *eos_table = ((pyEos*)object)->eos_table;
  return 1;
}

/* pyEos module */

PyObject* _loadeos(PyObject *self, PyObject *args) {
  char *datafile, *speciesfile;
  t_eos_table* eos_table;
  int oldformat, ret;
  
  oldformat = 0;
  if (!PyArg_ParseTuple( args, "ss|i:loadeos( datafile, speciesfile, [oldformat] )", &datafile, &speciesfile, &oldformat )) {
    return 0;
  }

  eos_table = (t_eos_table*)malloc( sizeof( t_eos_table ) );
        if (eos_table == NULL) {
          PyErr_SetFromErrno(PyExc_MemoryError);
          return NULL;
        }
  if (!oldformat) {
    ret = eos_init( eos_table, datafile, speciesfile );
  } else {
    ret = eos_init_old( eos_table, datafile, speciesfile );
  }

        if (ret != 0) {
          free(eos_table);
          PyErr_SetString(PyExc_RuntimeError, "error in EoS initialization");
          return NULL;
        }

  return PyGetEosObject( eos_table );
}

static PyMethodDef pyeosMethods[] = {
  { "loadeos", _loadeos, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "pyeos", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  pyeosMethods,        /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_pyeos(void)
{
  import_array();
  pyEosType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&pyEosType) < 0)
    return NULL;
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initpyeos(void)
{
  Py_InitModule( "pyeos", pyeosMethods );
  import_array();
}
#endif
