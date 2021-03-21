#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gadgetSnap.h"

#define G 6.672e-8 /* gravitational constant */

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

inline void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
  PyDict_SetItemString(dict, key, object);
  Py_DECREF(object);
}

PyArrayObject* createPyDoubleArray( double *data, int length ) {
  PyArrayObject* pyData;
  
  npy_intp dims[1];
  dims[0] = length;
  
  pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_DOUBLE );
  memcpy( PyArray_DATA(pyData), data, length*sizeof(double) );
  
  return pyData;
}

PyArrayObject* createPyIntArray( int *data, int length ) {
  PyArrayObject* pyData;
  
  npy_intp dims[1];
  dims[0] = length;
    
  pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_INT );
  memcpy( PyArray_DATA(pyData), data, length*sizeof(int) );
  
  return pyData;
}

PyObject* _getBoundMass(PyObject *self, PyObject *args) {
  char *filename;
  t_gadgetSnap *snap;
  int npart, i, j;
  PyObject *dict;
  
  double cx, cy, cz;
  double vx, vy, vz;
  int ncells, nstar, getbids;
  double radius, dr;
  double *massrad;
  
  float *pos, *vel, *mass, *pot;
  int *id;
  double r, v2;
  int cell, nbound;
  double bmass, bmass_star, bmass_sn, bmass2, bmass2_star, bmass2_sn;
  double bpos[3], bvel[3];
  int *bid_sn, *bid_sn_tmp, bid_sn_count, bid_sn_size;
  int has_pot;
  
  getbids = 0;
  ncells = 10000;
  if (!PyArg_ParseTuple( args, "siddddddd|ii:getBoundMass( snap, nstar, cx, cy, cz, vx, vy, vz, radius, [getbids, ncells]  )", &filename, &nstar, &cx, &cy, &cz, &vx, &vy, &vz, &radius, &getbids, &ncells )) {
    return 0;
  }
  printf( "vel: %g %g %g\n", vx, vy, vz );

  massrad = (double*)malloc( ncells * sizeof( double ) );
  memset( massrad, 0, ncells * sizeof(double) );
  dr = radius / ncells;
  
  snap = gc_init( filename );
  has_pot = gc_hasprop( snap, "POT " );

  gc_enableprop( snap, "POS " );
  gc_enableprop( snap, "MASS" );

  bid_sn = 0;
  
  /* get the radial mass profile first */
  gc_reset( snap );
  npart = gc_readnext( snap );
  while (npart > 0) {
    pos  = (float*)gc_getdata( snap, "POS " );
    mass = (float*)gc_getdata( snap, "MASS" );
    
    for (i=0; i<npart; i++) {
      r = sqrt( (pos[0]-cx)*(pos[0]-cx) + (pos[1]-cy)*(pos[1]-cy) + (pos[2]-cz)*(pos[2]-cz) );
      
      if (r < radius) {
        cell = floor( r / dr );
        massrad[ cell ] += mass[0];
      }
      
      pos += 3;
      mass += 1;
    }
    
    npart = gc_readnext( snap );
  }
  
  /* sum up masses, so that massrad is the total mass within a given radius */
  for (i=1; i<ncells; i++) {
    massrad[i] += massrad[i-1];
  }
  
  gc_enableprop( snap, "VEL " );
  gc_enableprop( snap, "ID  " );
  
  bmass = 0.0;
  bmass_star = 0.0;
  bmass_sn = 0.0;
  nbound = 0;
  
  for (j=0; j<3; j++) {
    bpos[j] = 0.0;
    bvel[j] = 0.0;
  }
  
  bid_sn = 0;
  if (getbids) {
    bid_sn_size = 1000;
    bid_sn = (int*)malloc( bid_sn_size * sizeof(double) );
    memset( bid_sn, 0, bid_sn_size * sizeof(double) );
    bid_sn_count = 0;
  }
  
  /* check which particles are bound to the star */
  gc_reset( snap );
  npart = gc_readnext( snap );
  while (npart > 0) {
    pos  = (float*)gc_getdata( snap, "POS " );
    vel  = (float*)gc_getdata( snap, "VEL " );
    mass = (float*)gc_getdata( snap, "MASS" );
    id   = (int*)gc_getdata( snap, "ID  " );
    
    for (i=0; i<npart; i++) {
      r = sqrt( (pos[0]-cx)*(pos[0]-cx) + (pos[1]-cy)*(pos[1]-cy) + (pos[2]-cz)*(pos[2]-cz) );
      
      if (r < radius) {
        v2 = (vel[0]-vx)*(vel[0]-vx) + (vel[1]-vy)*(vel[1]-vy) + (vel[2]-vz)*(vel[2]-vz);
        cell = floor( r / dr );
        
        if (v2 < 2. * G * massrad[cell] / r) {
          bmass += mass[0];
          if (id[0] <= nstar) {
            bmass_star += mass[0];
          } else {
            bmass_sn += mass[0];
            if (getbids) {
              if (bid_sn_count == bid_sn_size) {
                bid_sn_tmp = (int*)malloc( bid_sn_size * 2 * sizeof(int) );
                memcpy( bid_sn_tmp, bid_sn, bid_sn_size * sizeof(int) );
                free( bid_sn );
                bid_sn = bid_sn_tmp;
                bid_sn_size *= 2;
              }
              bid_sn[bid_sn_count] = id[0];
              bid_sn_count++;
            }
          }
          for (j=0; j<3; j++) {
            bpos[j] += pos[j];
            bvel[j] += vel[j];
          }
          nbound++;
        }
      }
      
      pos += 3;
      vel += 3;
      mass += 1;
      id += 1;
    }
    npart = gc_readnext( snap );
  }
  
  free( massrad );
  
  for (j=0; j<3; j++) {
    bpos[j] /= nbound;
    bvel[j] /= nbound;
  }

  if (has_pot) {
    for (j=0; j<3; j++) {
      bpos[j] = 0.0;
      bvel[j] = 0.0;
    }

    if (getbids) {
      if (bid_sn) free( bid_sn );
      bid_sn_size = 1000;
      bid_sn = (int*)malloc( bid_sn_size * sizeof(double) );
      memset( bid_sn, 0, bid_sn_size * sizeof(double) );
      bid_sn_count = 0;
    }

    gc_enableprop( snap, "POT " );
    bmass2 = 0.0;
    bmass2_star = 0.0;
    bmass2_sn = 0.0;
    nbound = 0;

    gc_reset( snap );
    npart = gc_readnext( snap );
    while (npart > 0) {
      pos  = (float*)gc_getdata( snap, "POS " );
      vel  = (float*)gc_getdata( snap, "VEL " );
      mass = (float*)gc_getdata( snap, "MASS" );
      id   = (int*)gc_getdata( snap, "ID  " );
      pot  = (float*)gc_getdata( snap, "POT " );
    
      for (i=0; i<npart; i++) {
        v2 = (vel[0]-vx)*(vel[0]-vx) + (vel[1]-vy)*(vel[1]-vy) + (vel[2]-vz)*(vel[2]-vz);
        
        if (v2 < -2.*pot[0]) {
                bmass2 += mass[0];
          if (id[0] <= nstar) {
            bmass2_star += mass[0];
          } else {
            bmass2_sn += mass[0];
            if (getbids) {
              if (bid_sn_count == bid_sn_size) {
                bid_sn_tmp = (int*)malloc( bid_sn_size * 2 * sizeof(int) );
                memcpy( bid_sn_tmp, bid_sn, bid_sn_size * sizeof(int) );
                free( bid_sn );
                bid_sn = bid_sn_tmp;
                bid_sn_size *= 2;
              }
              bid_sn[bid_sn_count] = id[0];
              bid_sn_count++;
            }
          }

          for (j=0; j<3; j++) {
            bpos[j] += pos[j];
            bvel[j] += vel[j];
          }
          nbound++;
        }
      
        pos += 3;
        vel += 3;
        mass += 1;
        id += 1;
        pot += 1;
      }
      npart = gc_readnext( snap );
    }

    for (j=0; j<3; j++) {
      bpos[j] /= nbound;
      bvel[j] /= nbound;
    }
  }

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "time", (PyObject*)PyFloat_FromDouble( snap->time ) );
  PyDict_SetStolenItem( dict, "bmass", (PyObject*)PyFloat_FromDouble( bmass ) );
  PyDict_SetStolenItem( dict, "bstar", (PyObject*)PyFloat_FromDouble( bmass_star ) );
  PyDict_SetStolenItem( dict, "bsn", (PyObject*)PyFloat_FromDouble( bmass_sn ) );
  if (has_pot) {
    PyDict_SetStolenItem( dict, "bmass2", (PyObject*)PyFloat_FromDouble( bmass2 ) );
    PyDict_SetStolenItem( dict, "bstar2", (PyObject*)PyFloat_FromDouble( bmass2_star ) );
    PyDict_SetStolenItem( dict, "bsn2", (PyObject*)PyFloat_FromDouble( bmass2_sn ) );
  }
  PyDict_SetStolenItem( dict, "bpos", (PyObject*)createPyDoubleArray( bpos, 3 ) );
  PyDict_SetStolenItem( dict, "bvel", (PyObject*)createPyDoubleArray( bvel, 3 ) );

  if (getbids) {
    PyDict_SetStolenItem( dict, "bid_sn", (PyObject*)createPyIntArray( bid_sn, bid_sn_count ) );
    free( bid_sn );
  }
  
  gc_deinit( snap );
  
  return dict;
}

static PyMethodDef boundmassmethods[] = {
  { "getboundmass", _getBoundMass, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "boundmass", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  boundmassmethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_boundmass(void)
{
  import_array();
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initboundmass(void)
{
  Py_InitModule( "boundmass", boundmassmethods );
  import_array();
}
#endif
