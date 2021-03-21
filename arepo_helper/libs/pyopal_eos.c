#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pyopal_eos.h"
#include "opal_eos.h"

typedef struct {
	PyObject_HEAD;
	struct opal_eos_table *opal_eos_table;
} pyOpalEos;

/* forward declarations */
void pyOpalEosDealloc( PyObject* self );
static PyObject* pyOpalEos_str( pyOpalEos* self );

void pyOpalEosDealloc( PyObject* self ) {
	opaleos_deinit( ((pyOpalEos*)self)->opal_eos_table );
	PyObject_Del( self );
}

PyObject* pyOpalEos_str( pyOpalEos* self ) {
	return PyUnicode_FromString( "opal_eos" );
}

PyObject* pyOpalEos_getxlim( pyOpalEos* self, PyObject* args ) {
	return Py_BuildValue( "{s:d,s:d}", "minx", self->opal_eos_table->minx, "maxx", self->opal_eos_table->maxx);
}

PyObject* pyOpalEos_getrholim( pyOpalEos* self, PyObject* args ) {
	return Py_BuildValue( "{s:d,s:d}", "minrho", self->opal_eos_table->minrho, "maxrho", self->opal_eos_table->maxrho);
}

PyObject* pyOpalEos_getrhovalues( pyOpalEos* self, PyObject* args ) {
        int i;
        PyListObject * rhos = (PyListObject *)PyList_New(NRHO);
        for (i = 0; i < NRHO; i++)
          PyList_SET_ITEM(rhos, i, PyFloat_FromDouble(self->opal_eos_table->rho[i]));
        
	return Py_BuildValue( "O", rhos);
}

PyObject* pyOpalEos_gettvalues( pyOpalEos* self, PyObject* args ) {
        int i;
        PyListObject * ts = (PyListObject *)PyList_New(NT);
        for (i = 0; i < NT; i++)
          PyList_SET_ITEM(ts, i, PyFloat_FromDouble(self->opal_eos_table->t[i]*1e6));
        
	return Py_BuildValue( "O", ts);
}

PyObject* pyOpalEos_gettlim( pyOpalEos* self, PyObject* args ) {
        int i;
        PyListObject * mints = (PyListObject *)PyList_New(NRHO);
        for (i = 0; i < NRHO; i++)
          PyList_SET_ITEM(mints, i, PyFloat_FromDouble(self->opal_eos_table->mint[i]*1e6));
        
	return Py_BuildValue( "{s:O,s:d}", "mint", mints, "maxt", self->opal_eos_table->maxt*1e6);
}

PyObject* pyOpalEos_tgiven( pyOpalEos* self, PyObject* args ) {
	double x, rho, temp;
	struct opal_eos_result res;

	if ( !PyArg_ParseTuple( args, "ddd:opaleos.tgiven(x, rho, temp)", &x, &rho, &temp ) ) {
		return 0;
	}
	
	if (opaleos_all_tgiven( self->opal_eos_table, x, rho, temp, &res ) < 0)
          {
              //PyErr_SetString( PyExc_ValueError, "Error in OPAL EOS.\n" );
              //return 0;
	  }
	return Py_BuildValue( "{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}", "e", res.e, "p", res.p, 
            "s", res.s, "dedrho", res.dedrho, "cv", res.cv, "chit", res.chit, "chir", res.chir, 
            "gamma1", res.gamma1, "gamma2", res.gamma2, "gamma3", res.gamma3 );
}

PyObject* pyOpalEos_egiven( pyOpalEos* self, PyObject* args ) {
	double x, rho, temp, e;
	struct opal_eos_result res;

	if ( !PyArg_ParseTuple( args, "ddd:opaleos.egiven(x, rho, e)", &x, &rho, &e ) ) {
		return 0;
	}
	
	
        temp = -1.0; /* let egiven guess T */
	if (opaleos_all_egiven( self->opal_eos_table, x, rho, e, &temp, &res ) < 0)
          {
              //PyErr_SetString( PyExc_ValueError, "Error in OPAL EOS.\n" );
              //return 0;
	  }
	return Py_BuildValue( "{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}", "t", temp, "p", res.p, 
            "s", res.s, "dedrho", res.dedrho, "cv", res.cv, "chit", res.chit, "chir", res.chir, 
            "gamma1", res.gamma1, "gamma2", res.gamma2, "gamma3", res.gamma3 );
}

PyObject* pyOpalEos_pgiven( pyOpalEos* self, PyObject* args ) {
	double x, rho, temp, p;
	struct opal_eos_result res;

	if ( !PyArg_ParseTuple( args, "ddd:opaleos.pgiven(x, rho, p)", &x, &rho, &p ) ) {
		return 0;
	}
	
	
        temp = -1.0; /* let pgiven guess T */
	if (opaleos_all_pgiven( self->opal_eos_table, x, rho, p, &temp, &res ) < 0)
          {
              //PyErr_SetString( PyExc_ValueError, "Error in OPAL EOS.\n" );
              //return 0;
	  }
	return Py_BuildValue( "{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}", "t", temp, "e", res.e, 
            "s", res.s, "dedrho", res.dedrho, "cv", res.cv, "chit", res.chit, "chir", res.chir, 
            "gamma1", res.gamma1, "gamma2", res.gamma2, "gamma3", res.gamma3, "p", res.p );
}

PyObject* pyOpalEos_ptgiven( pyOpalEos* self, PyObject* args ) {
	double x, temp, p;
	struct opal_eos_result res;

	if ( !PyArg_ParseTuple( args, "ddd:opaleos.pgiven(x, p, t)", &x, &p, &temp ) ) {
		return 0;
	}
	
	
	if (opaleos_all_ptgiven( self->opal_eos_table, x, p, temp, &res ) < 0)
          {
              //PyErr_SetString( PyExc_ValueError, "Error in OPAL EOS.\n" );
              //return 0;
	  }
	return Py_BuildValue( "{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}", "rho", res.rho, "e", res.e, 
            "s", res.s, "dedrho", res.dedrho, "cv", res.cv, "chit", res.chit, "chir", res.chir, 
            "gamma1", res.gamma1, "gamma2", res.gamma2, "gamma3", res.gamma3, "p", res.p );
}


static PyMethodDef pyOpalEosMethods[] =
{
{ "egiven", (getattrofunc)pyOpalEos_egiven, METH_VARARGS, 
  "egiven(x, rho, e)\nReturns dictionary with all EOS quantities for given x, rho, e.\n(x: mass fraction of H; rho: density; e: internal energy)" },
{ "pgiven", (getattrofunc)pyOpalEos_pgiven, METH_VARARGS, 
  "pgiven(x, rho, p)\nReturns dictionary with all EOS quantities for given x, rho, p.\n(x: mass fraction of H; rho: density; p: pressure)" },
{ "tgiven", (getattrofunc)pyOpalEos_tgiven, METH_VARARGS, 
  "tgiven(x, rho, t)\nReturns dictionary with all EOS quantities for given x, rho, t.\n(x: mass fraction of H; rho: density; t: temperature)" },
{ "ptgiven", (getattrofunc)pyOpalEos_ptgiven, METH_VARARGS, 
  "ptgiven(x, p, t)\nReturns dictionary with all EOS quantities for given x, p, t.\n(x: mass fraction of H; p: pressure; t: temperature)" },
{ "getxlim", (getattrofunc)pyOpalEos_getxlim, METH_VARARGS, 
  "getxlim()\nReturns limits of table in X direction." },
{ "getrholim", (getattrofunc)pyOpalEos_getrholim, METH_VARARGS, 
  "getrholim()\nReturns limits of table in rho direction." },
{ "gettlim", (getattrofunc)pyOpalEos_gettlim, METH_VARARGS, 
  "gettlim()\nReturns limits of table in T direction." },
{ "getrhovalues", (getattrofunc)pyOpalEos_getrhovalues, METH_VARARGS, 
  "getrhovalues()\nReturns values of rho for table." },
{ "gettvalues", (getattrofunc)pyOpalEos_gettvalues, METH_VARARGS, 
  "gettvalues()\nReturns values of T for table." },
{ NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3
static PyTypeObject pyOpalEosType =
{
PyVarObject_HEAD_INIT( NULL, 0)
"pyopaleos",
sizeof( pyOpalEos ),
0,
pyOpalEosDealloc,
0,
0,
0,0,0,0,0,0,0,0,
( reprfunc ) pyOpalEos_str,
0,0,0,0,0,0,0,0,0,0,
0,
pyOpalEosMethods, /* tp_methods */
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
};
#else
static PyObject* pyOpalEos_getAttr( pyOpalEos* self, char* name );

static PyTypeObject pyOpalEosType =
{
PyObject_HEAD_INIT( &PyType_Type )
0,
"pyopaleos",
sizeof( pyOpalEosType ),
0,
pyOpalEosDealloc,
0,
( getattrfunc ) pyOpalEos_getAttr,
0,0,0,0,0,0,0,0,
( reprfunc ) pyOpalEos_str,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
};

PyObject* pyOpalEos_getAttr( pyOpalEos* self, char* name ) {
  return Py_FindMethod( pyOpalEosMethods, (PyObject*)self, name );
}
#endif


PyObject* PyGetOpalEosObject( struct opal_eos_table* opal_eos_table ) {
	pyOpalEos* eos;
	
	if (opal_eos_table) {
		eos = (pyOpalEos*)PyObject_New( pyOpalEos, &pyOpalEosType );
		eos->opal_eos_table = opal_eos_table;
		return (PyObject*)eos;
	} else {
		Py_RETURN_NONE;
	}
}

int pyConvertOpalEos( PyObject* object, struct opal_eos_table** opal_eos_table ) {
	if (strcmp(object->ob_type->tp_name, "pyopaleos" )) {
		PyErr_BadArgument();
		return 0;
	}
	
	*opal_eos_table = ((pyOpalEos*)object)->opal_eos_table;
	return 1;
}

/* pyEos module */

PyObject* _loadopaleos(PyObject *self, PyObject *args) {
	char *datafile;
	struct opal_eos_table* opal_eos_table;
	
	if (!PyArg_ParseTuple( args, "s:loadopal_eos( datafile )", &datafile )) {
		return 0;
	}

	opal_eos_table = (struct opal_eos_table*)malloc( sizeof( struct opal_eos_table ) );
	opal_eos_table = opaleos_init( datafile );
	
	return PyGetOpalEosObject( opal_eos_table );
}

static PyMethodDef pyopaleosMethods[] = {
	{ "loadopal_eos", _loadopaleos, METH_VARARGS, 
          "loadopal_eos(datafile)\nLoads OPAL EOS from specified data file (usually EOS5_data)." },
	{ NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "pyopal_eos", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  pyopaleosMethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_pyopal_eos(void)
{
	import_array();
        pyOpalEosType.tp_new = PyType_GenericNew;
        if (PyType_Ready(&pyOpalEosType) < 0)
          return NULL;
	return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initpyopal_eos(void)
{
	Py_InitModule( "pyopal_eos", pyopaleosMethods );
	import_array();
}
#endif
