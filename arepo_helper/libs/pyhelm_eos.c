#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pyhelm_eos.h"
#include "helm_eos.h"

static void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
  PyDict_SetItemString(dict, key, object);
  Py_DECREF(object);
}

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
{ NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3
static PyTypeObject pyHelmEosType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "pyhelmeos",               /* tp_name */
    sizeof(pyHelmEos),         /* tp_basicsize */
    0,                         /* tp_itemsize */
    pyHelmEosDealloc,                         /* tp_dealloc */
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
    (reprfunc)pyHelmEos_str,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Helmholz EOS",            /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    pyHelmEosMethods,          /* tp_methods */
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
static PyObject* pyHelmEos_getAttr( pyHelmEos* self, char* name );

static PyTypeObject pyHelmEosType =
{
PyObject_HEAD_INIT( &PyType_Type )
0,
"pyhelmeos",
sizeof( pyHelmEosType ),
0,
pyHelmEosDealloc,
0,
( getattrfunc ) pyHelmEos_getAttr,
0,0,0,0,0,0,0,0,
( reprfunc ) pyHelmEos_str,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
};

PyObject* pyHelmEos_getAttr( pyHelmEos* self, char* name ) {
  return Py_FindMethod( pyHelmEosMethods, (PyObject*)self, name );
}
#endif

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
  struct eos_result res;

  if ( !PyArg_ParseTuple( args, "dO!d:helmeos.egiven(rho, xnuc, e)", &rho, &PyArray_Type, &pyXnuc, &e ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
    printf( "%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies );
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
  }
  
  temp = -1.;
  eos_calc_egiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), e, &temp, &res );
  
  return Py_BuildValue( "dd", temp, res.p.v );
}

PyObject* pyHelmEos_pgiven( pyHelmEos* self, PyObject* args ) {
  PyArrayObject *pyXnuc;
  double rho, p, temp;
  struct eos_result res;

  if ( !PyArg_ParseTuple( args, "dO!d:helmeos.pgiven(rho, xnuc, p)", &rho, &PyArray_Type, &pyXnuc, &p ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
    printf( "%d %d\n", PyArray_NDIM(pyXnuc), self->helm_eos_table->nspecies );
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
  }
  
  temp = -1.;
  eos_calc_pgiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), p, &temp, &res );
  
  return Py_BuildValue( "dd", temp, res.e.v );
}

PyObject* pyHelmEos_tgiven( pyHelmEos* self, PyObject* args ) {
  PyArrayObject *pyXnuc;
  double rho, temp;
  struct eos_result res;

  if ( !PyArg_ParseTuple( args, "dO!d:helmeos.tgiven(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
  }
  
  eos_calc_tgiven( self->helm_eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &res );
  
  return Py_BuildValue( "dddd", res.e.v, res.e.dtemp, res.p.v, res.sound );
}

PyObject* pyHelmEos_tgivenfull( pyHelmEos* self, PyObject* args ) {
  PyArrayObject *pyXnuc;
  double rho, temp;
  struct eos_result res;

  if ( !PyArg_ParseTuple( args, "dO!d:helmeos.tgivenfull(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp ) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
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
  struct eos_result res;

        rho = -1;
  if ( !PyArg_ParseTuple( args, "dO!d|d:helmeos.ptgivenfull(p, xnuc, temp, [rho])", &p, &PyArray_Type, &pyXnuc, &temp, &rho) ) {
    return 0;
  }
  
  if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->helm_eos_table->nspecies) {
    PyErr_SetString( PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile." );
    return 0;
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
  if (strcmp(object->ob_type->tp_name, "pyhelmeos" )) {
    PyErr_BadArgument();
    return 0;
  }
  
  *helm_eos_table = ((pyHelmEos*)object)->helm_eos_table;
  return 1;
}

/* pyEos module */

PyObject* _eos(PyObject *self, PyObject *args) {
  char *datafile, *speciesfile;
  t_helm_eos_table* helm_eos_table;
  
  int useCoulombCorrections = 0;
  int verbose = 0;
  if (!PyArg_ParseTuple( args, "ss|ii:loadhelmeos( datafile, speciesfile, [useCoulombCorrections,verbose] )", &datafile, &speciesfile, &useCoulombCorrections, &verbose )) {
    return 0;
  }

  helm_eos_table = (t_helm_eos_table*)malloc( sizeof( t_helm_eos_table ) );
  eos_init( helm_eos_table, datafile, speciesfile, useCoulombCorrections, verbose );
  
  return PyGetHelmEosObject( helm_eos_table );
}

/* pyEos module */

PyObject* _loadhelmeos(PyObject *self, PyObject *args) {
  char *datafile, *speciesfile;
  t_helm_eos_table* helm_eos_table;

  int useCoulombCorrections = 0;
  int verbose = 0;
  if (!PyArg_ParseTuple( args, "ss|ii:loadhelmeos( datafile, speciesfile, [useCoulombCorrections,verbose] )", &datafile, &speciesfile, &useCoulombCorrections, &verbose )) {
    return 0;
  }

  helm_eos_table = (t_helm_eos_table*)malloc( sizeof( t_helm_eos_table ) );
  eos_init( helm_eos_table, datafile, speciesfile, useCoulombCorrections, verbose );

  return PyGetHelmEosObject( helm_eos_table );
}


static PyMethodDef pyhelmeosMethods[] = {
  { "loadhelm_eos", _loadhelmeos, METH_VARARGS, 
    "loadhelm_eos(datafile, speciesfile, [useCoulombCorrections, verbose]).\nLoad Helmholz EOS from specified datafile and speciesfile." },
  { NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "pyhelm_eos", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  pyhelmeosMethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_pyhelm_eos(void)
{
  import_array();
  pyHelmEosType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&pyHelmEosType) < 0)
    return NULL;
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initpyhelm_eos(void)
{
  Py_InitModule( "pyhelm_eos", pyhelmeosMethods );
  import_array();
}
#endif
