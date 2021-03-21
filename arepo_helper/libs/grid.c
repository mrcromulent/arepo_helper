#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define sqr(x) ((x)*(x))


inline void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
  PyDict_SetItemString(dict, key, object);
  Py_DECREF(object);
}

PyArrayObject* createPyArray( double *data, int dim1, int dim2 ) {
  PyArrayObject* pyData;

  if (dim1 == 1 || dim2 == 1) {
    npy_intp dims[1];
    dims[0] = dim1 * dim2;
    pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_DOUBLE );
    memcpy( PyArray_DATA(pyData), data, dim1*dim2*sizeof(double) );
  } else {
    npy_intp dims[2];
    dims[0] = dim1;
    dims[1] = dim2;
    pyData = (PyArrayObject *)PyArray_SimpleNew( 2, dims, NPY_DOUBLE );
    memcpy( PyArray_DATA(pyData), data, dim1*dim2*sizeof(double) );
  }   
  
  return pyData;
}

PyObject* _calcdmdv(PyObject *self, PyObject *args) {
  PyObject *dict;
  PyArrayObject *pyMass, *pyVel, *pyX, *pyY, *pyZ;
  int ncells;
  double dv, dmax, massloss;
  double vel, *mass, *v;
  float *x, *y, *z;
  int nx, ny, nz, i, j, k, cell, dir;

  dv   =   1e8;
  dmax = 1.5e9;
  dir = 0;
  if (!PyArg_ParseTuple( args, "O!O!O!O!O!|ddi:calcdmdv( mass, vel, x, y, z, [dv,dmax,dir]  )", &PyArray_Type, &pyMass, &PyArray_Type, &pyVel, &PyArray_Type, &pyX, &PyArray_Type, &pyY, &PyArray_Type, &pyZ, &dv, &dmax, &dir )) {
    return 0;
  }

  if (PyArray_NDIM(pyMass) != 3 || PyArray_TYPE(pyMass) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "rho has to be of dimensions [nx,ny,nz] and type double" );
    return 0;
  }

  nx = PyArray_DIMS(pyMass)[0];
  ny = PyArray_DIMS(pyMass)[1];
  nz = PyArray_DIMS(pyMass)[2];

  if (PyArray_NDIM(pyVel) != 3 || PyArray_TYPE(pyVel) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "vel has to be of dimensions [nx,ny,nz] and type double" );
    return 0;
  }

  if (dir < -1 || dir > 1) {
    PyErr_SetString( PyExc_ValueError, "if set, dir has to bei -1, 0 or 1" );
    return 0;
  }

  if (PyArray_NDIM(pyX) != 1 || PyArray_TYPE(pyX) != NPY_FLOAT || PyArray_DIMS(pyX)[0] != nx) {
    PyErr_SetString( PyExc_ValueError, "x has to be of dimensions [nx] and type float" );
    return 0;
  }

  if (PyArray_NDIM(pyY) != 1 || PyArray_TYPE(pyY) != NPY_FLOAT || PyArray_DIMS(pyY)[0] != ny) {
    PyErr_SetString( PyExc_ValueError, "y has to be of dimensions [ny] and type float" );
    return 0;
  }

  if (PyArray_NDIM(pyZ) != 1 || PyArray_TYPE(pyZ) != NPY_FLOAT || PyArray_DIMS(pyZ)[0] != nz) {
    PyErr_SetString( PyExc_ValueError, "z has to be of dimensions [nz] and type float" );
    return 0;
  }

  x = (float*)PyArray_DATA(pyX);
  y = (float*)PyArray_DATA(pyY);
  z = (float*)PyArray_DATA(pyZ);

  ncells = ceil( dmax / dv );
  mass = (double*)malloc( ncells * sizeof(double) );
  memset( mass, 0, ncells * sizeof(double) );

  massloss = 0.0;

  for (i=0; i<nx; i++) for (j=0; j<ny; j++) for (k=0; k<nz; k++) {
    vel = *(double*)( PyArray_DATA(pyVel) + PyArray_STRIDES(pyVel)[0]*i + PyArray_STRIDES(pyVel)[1]*j + PyArray_STRIDES(pyVel)[2]*k );
    cell = floor( vel / dv + 0.5 );

    if (cell < ncells) {
        if ((dir == 0) || (dir == -1 && z[k]*z[k] <= x[i]*x[i] + y[j]*y[j]) || (dir == 1 && z[k]*z[k] > x[i]*x[i] + y[j]*y[j]))  {
        mass[cell] += *(double*)( PyArray_DATA(pyMass) + PyArray_STRIDES(pyMass)[0]*i + PyArray_STRIDES(pyMass)[1]*j + PyArray_STRIDES(pyMass)[2]*k );
      }
    } else {
      massloss += *(double*)( PyArray_DATA(pyMass) + PyArray_STRIDES(pyMass)[0]*i + PyArray_STRIDES(pyMass)[1]*j + PyArray_STRIDES(pyMass)[2]*k );
    }
  }
  
  v = (double*)malloc( ncells * sizeof(double) );
  for (i=0; i<ncells; i++) v[i] = dv * i;

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "m", (PyObject*)createPyArray( mass, ncells, 1 ) );
  free( mass );
  PyDict_SetStolenItem( dict, "v", (PyObject*)createPyArray( v, ncells, 1 ) );
  free( v );

  return dict;
}

PyObject* _calcabund(PyObject *self, PyObject *args) {
  PyObject *result;
  PyArrayObject *pyDens, *pyRho, *pyAbund;
  int ncells, nrho, nspecies;
  double *dens, *rho, *sum;
  double rhomin, rhomax, rr;
  int i, j, cell;
  
  if (!PyArg_ParseTuple( args, "O!O!O!:calcabund( dens, rho, abund )", &PyArray_Type, &pyDens, &PyArray_Type, &pyRho, &PyArray_Type, &pyAbund )) {
    return 0;
  }

  if (PyArray_NDIM(pyDens) != 3 || PyArray_TYPE(pyDens) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "dens has to be of dimensions [nx,ny,nz] and type double" );
    return 0;
  }
  ncells = PyArray_DIMS(pyDens)[0] * PyArray_DIMS(pyDens)[1] * PyArray_DIMS(pyDens)[2];
  dens = (double*)PyArray_DATA(pyDens);

  if (PyArray_NDIM(pyRho) != 1 || PyArray_TYPE(pyDens) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "rho has to be of dimensions [nx] and type double" );
    return 0;
  }
  nrho = PyArray_DIMS(pyRho)[0];
  rho = (double*)PyArray_DATA(pyRho);

  if (PyArray_NDIM(pyAbund) != 2 || PyArray_TYPE(pyDens) != NPY_DOUBLE ) {
    PyErr_SetString( PyExc_ValueError, "abund has to be of dimensions [nx,ny] and type double" );
    return 0;
  }
  nspecies = PyArray_DIMS(pyAbund)[1];

  if (PyArray_DIMS(pyAbund)[0] != nrho) {
    PyErr_SetString( PyExc_ValueError, "abund and rho have to have the same amount of entries in the first dimension" );
    return 0;
  }

  rhomin = 1e20;
  for (i=0; i<nrho; i++) if (rho[i] < rhomin) rhomin = rho[i];
  rhomax = 0;
  for (i=0; i<nrho; i++) if (rho[i] > rhomax) rhomax = rho[i];

  sum = (double*)malloc( nspecies * sizeof(double) );
  memset( sum, 0, nspecies * sizeof(double) );

  for (i=0; i<ncells; i++) {
    rr = dens[i];
    
    if (rr <= rhomin || rr >= rhomax) continue;
    
    cell = 0;
    while (rho[cell] > rr)
      cell++;

    if (cell > 0 && abs(rho[cell]-rr) > abs(rho[cell-1]-rr)) cell--;

    for (j=0; j<nspecies; j++) {
      sum[j] += *(double*)( PyArray_DATA(pyAbund) + PyArray_STRIDES(pyAbund)[0]*cell + PyArray_STRIDES(pyAbund)[1]*j ) * rr;
    }
  }

  result = (PyObject*)createPyArray( sum, nspecies, 1 );
  free( sum );
  return result;
}

PyObject* _calcdrhodv(PyObject *self, PyObject *args) {
  PyObject *dict;
  PyArrayObject *pyRho, *pyVel;
  int ncells;
  double dv, dmax;
  double *v, *rhosum, vel, dens;
  int *rhocount;
  int nx, ny, nz, i, j, k, cell;

  dv   =   1e8;
  dmax = 1.5e9;
  if (!PyArg_ParseTuple( args, "O!O!|dd:calcdrhodv( rho, vel, [dv,dmax]  )", &PyArray_Type, &pyRho, &PyArray_Type, &pyVel, &dv, &dmax )) {
    return 0;
  }

  if (PyArray_NDIM(pyRho) != 3 || PyArray_TYPE(pyRho) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "rho has to be of dimensions [nx,ny,nz] and type double" );
    return 0;
  }

  nx = PyArray_DIMS(pyRho)[0];
  ny = PyArray_DIMS(pyRho)[1];
  nz = PyArray_DIMS(pyRho)[2];

  if (PyArray_NDIM(pyVel) != 3 || PyArray_TYPE(pyVel) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "vel has to be of dimensions [nx,ny,nz] and type double" );
    return 0;
  }

  ncells = ceil( dmax / dv );
  rhosum = (double*)malloc( ncells * sizeof(double) );
  memset( rhosum, 0, ncells * sizeof(double) );
  rhocount = (int*)malloc( ncells * sizeof(int) );
  memset( rhocount, 0, ncells * sizeof(int) );

  for (i=0; i<nx; i++) for (j=0; j<ny; j++) for (k=0; k<nz; k++) {
    vel = *(double*)( PyArray_DATA(pyVel) + PyArray_STRIDES(pyVel)[0]*i + PyArray_STRIDES(pyVel)[1]*j + PyArray_STRIDES(pyVel)[2]*k );
    cell = floor( vel / dv + 0.5 );

    if (cell < ncells) {
      dens = *(double*)( PyArray_DATA(pyRho) + PyArray_STRIDES(pyRho)[0]*i + PyArray_STRIDES(pyRho)[1]*j + PyArray_STRIDES(pyRho)[2]*k );
      if (dens > 1e-4) {
        rhosum[cell] += dens;
        rhocount[cell]++;
      }
    }
  }
  
  v = (double*)malloc( ncells * sizeof(double) );
  for (i=0; i<ncells; i++) v[i] = dv * i;
  for (i=0; i<ncells; i++) rhosum[i] /= rhocount[i];
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( rhosum, ncells, 1 ) );
  free( rhosum );
  PyDict_SetStolenItem( dict, "v", (PyObject*)createPyArray( v, ncells, 1 ) );
  free( v );

  free( rhocount );
  return dict;
}

static PyMethodDef gridMethods[] = {
  { "calcdmdv", _calcdmdv, METH_VARARGS, "" },
  { "calcabund", _calcabund, METH_VARARGS, "" },
  { "calcdrhodv", _calcdrhodv, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "grid", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  gridMethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_grid(void)
{
  import_array();
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initgrid(void)
{
  Py_InitModule( "grid", gridMethods );
  import_array();
}
#endif
