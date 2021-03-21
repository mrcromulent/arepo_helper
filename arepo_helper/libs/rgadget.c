#include <Python.h>
#include <arrayobject.h>

#include "rgadget.h"
#include "gadgetSnap.h"

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

PyArrayObject* createPyIntArray( int *data, int length ) {
  PyArrayObject* pyData;
  
  npy_intp dims[1];
  dims[0] = length;
  
  pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_INT );
  memcpy( PyArray_DATA(pyData), data, length*sizeof(int) );
  
  return pyData;
}

void gadget_writeHeader( FILE* fd, int npartall, int npart, int num_files ) {
  int header, footer, blocklength;
  char *name;
  int nparts[6], npartsall[6];
  double masses[6];
  double time, redshift;
  int flag_sfr, flag_feedback, flag_cooling;
  char dummy[128];

  header = footer = 8;
  blocklength = 264; /* 256 (data) + 8 (header+footer) */
  name = "HEAD";

  fwrite( &header, sizeof(int), 1, fd );
  fwrite( name, sizeof(char), 4, fd );
  fwrite( &blocklength, sizeof(int), 1, fd );
  fwrite( &footer, sizeof(int), 1, fd );

  memset( nparts, 0, 6*sizeof(int) );
  memset( npartsall, 0, 6*sizeof(int) );
  memset( masses, 0, 6*sizeof(double) );
  memset( dummy, 0, 128 );

  header = footer = 256;
  nparts[0] = npart;
  npartsall[0] = npartall;
  time = 0;
  redshift = 0;
  flag_sfr = 0;
  flag_feedback = 0;
  flag_cooling = 0;
  
  fwrite( &header, sizeof(int), 1, fd );
  fwrite( nparts, sizeof(int), 6, fd );
  fwrite( masses, sizeof(double), 6, fd );
  fwrite( &time, sizeof(double), 1, fd );
  fwrite( &redshift, sizeof(double), 1, fd );
  fwrite( &flag_sfr, sizeof(int), 1, fd );
  fwrite( &flag_feedback, sizeof(int), 1, fd );
  fwrite( npartsall, sizeof(int), 6, fd );
  fwrite( &flag_cooling, sizeof(int), 1, fd );
  fwrite( &num_files, sizeof(int), 1, fd );
  fwrite( dummy, sizeof(char), 128, fd );
  fwrite( &footer, sizeof(int), 1, fd );
}

void gadget_writeBlockFloat( FILE* fd, char *name, int npart, int dim, float *data ) {
  int header, footer, blocklength;

  header = footer = 8;
  blocklength = npart * dim * sizeof(float) + 8;

  fwrite( &header, sizeof(int), 1, fd );
  fwrite( name, sizeof(char), 4, fd );
  fwrite( &blocklength, sizeof(int), 1, fd );
  fwrite( &footer, sizeof(int), 1, fd );

  header = footer = npart * dim * sizeof(float);
  fwrite( &header, sizeof(int), 1, fd );
  fwrite( data, sizeof(float), npart*dim, fd );
  fwrite( &footer, sizeof(int), 1, fd );  
}

void gadget_writeBlockInt( FILE* fd, char *name, int npart, int dim, int *data ) {
  int header, footer, blocklength;

  header = footer = 8;
  blocklength = npart * dim * sizeof(int) + 8;

  fwrite( &header, sizeof(int), 1, fd );
  fwrite( name, sizeof(char), 4, fd );
  fwrite( &blocklength, sizeof(int), 1, fd );
  fwrite( &footer, sizeof(int), 1, fd );

  header = footer = npart * dim * sizeof(float);
  fwrite( &header, sizeof(int), 1, fd );
  fwrite( data, sizeof(int), npart*dim, fd );
  fwrite( &footer, sizeof(int), 1, fd );
}

PyObject* _mergeICs(PyObject *self, PyObject *args) {
  FILE *fout;
  t_gadgetSnap *snap1, *snap2;
  char *fileout, *file1, *file2, fn[255];
  PyArrayObject *off1, *off2, *voff1, *voff2;
  double *data_off1, *data_off2, *data_voff1, *data_voff2, scale1, scale2;
  int npartalltot, start, i, j, split, nspecies, filenr, filenr1, filenr2, files;
  int *ids;
  float *value;

  off1 = off2 = 0;
  voff1 = voff2 = 0;
  scale1 = scale2 = 1.0;
  if (!PyArg_ParseTuple( args, "sss|O!O!O!O!dd:mergeICs( fileout, filename1, filename2, [offset1, offset2, voffset1, voffset2, scale1, scale2] )", &fileout, &file1, &file2, &PyArray_Type, &off1, &PyArray_Type, &off2, &PyArray_Type, &voff1, &PyArray_Type, &voff2, &scale1, &scale2 )) {
    Py_RETURN_FALSE;
  }

  snap1 = gc_init( file1 );
  snap2 = gc_init( file2 );

  if (!snap1 || !snap2) Py_RETURN_FALSE;

  if (gc_getdim( snap1, "XNUC" ) != gc_getdim( snap2, "XNUC" )) {
    printf( "Both snapshots have to have the same number of nuclear species.\n" );
    Py_RETURN_FALSE;
  } else {
    nspecies = gc_getdim( snap1, "XNUC" );
  }

  data_off1 = data_off2 = data_voff1 = data_voff2 = NULL;

  if (off1) data_off1 = (double*)PyArray_DATA(off1);
  if (off2) data_off2 = (double*)PyArray_DATA(off2);
  if (voff1) data_voff1 = (double*)PyArray_DATA(voff1);
  if (voff2) data_voff2 = (double*)PyArray_DATA(voff2);
  
  npartalltot = snap1->npartalltot + snap2->npartalltot;
  printf( "Total number of particles is %d (%d, %d)\n", npartalltot, snap1->npartalltot, snap2->npartalltot );
  
  if (snap1->multipleFiles > 1 || snap2->multipleFiles > 1 || npartalltot > MAXPARTPERICFILE) {
    split = 1;
  } else {
    split = 0;
  }

  if (!split) {
    fout = fopen( fileout, "w" );

    gadget_writeHeader( fout, npartalltot, npartalltot, 0 );
    snap1->particlesLoadAtOnce = snap1->npartalltot;
    snap2->particlesLoadAtOnce = snap2->npartalltot;

    /* POS and VEL */
    value = (float*)malloc( 3 * npartalltot * sizeof(float) );

    gc_readprop( snap1, "POS ", 0, snap1->nparttot );
    memcpy( value, gc_getdata( snap1 , "POS " ), 3 * snap1->nparttot * sizeof(float) );
    gc_freeprop( snap1, "POS " );

    gc_readprop( snap2, "POS ", 0, snap2->nparttot );
    memcpy( &value[3*snap1->nparttot], gc_getdata( snap2 , "POS " ), 3 * snap2->nparttot * sizeof(float) );
    gc_freeprop( snap2, "POS " );

    /* apply offsets */
    start = 3 * snap1->nparttot;
    if (off1) for (i=0; i<snap1->nparttot; i++) for (j=0; j<3; j++) value[i*3+j]       = value[i*3+j]*scale1 + data_off1[j];
    if (off2) for (i=0; i<snap2->nparttot; i++) for (j=0; j<3; j++) value[start+i*3+j] = value[start+i*3+j]*scale2 + data_off2[j];

    gadget_writeBlockFloat( fout, "POS ", npartalltot, 3, value );

    gc_readprop( snap1, "VEL ", 0, snap1->nparttot );
    memcpy( value, gc_getdata( snap1 , "VEL " ), 3 * snap1->nparttot * sizeof(float) );
    gc_freeprop( snap1, "VEL " );

    gc_readprop( snap2, "VEL ", 0, snap2->nparttot );
    memcpy( &value[3*snap1->nparttot], gc_getdata( snap2 , "VEL " ), 3 * snap2->nparttot * sizeof(float) );
    gc_freeprop( snap2, "VEL " );

    /* apply offsets */
    start = 3 * snap1->nparttot;
    if (voff1) for (i=0; i<snap1->nparttot; i++) for (j=0; j<3; j++) value[i*3+j] += data_voff1[j];
    if (voff2) for (i=0; i<snap2->nparttot; i++) for (j=0; j<3; j++) value[start+i*3+j] += data_voff2[j];

    gadget_writeBlockFloat( fout, "VEL ", npartalltot, 3, value );

    free( value );

    /* IDs */
    ids = (int*)malloc( npartalltot * sizeof(int) );

    gc_readprop( snap1, "ID  ", 0, snap1->nparttot );
    memcpy( ids, gc_getdata( snap1 , "ID  " ), snap1->nparttot * sizeof(int) );
    gc_freeprop( snap1, "ID  " );

    gc_readprop( snap2, "ID  ", 0, snap2->nparttot );
    memcpy( &ids[snap1->nparttot], gc_getdata( snap2 , "ID  " ), snap2->nparttot * sizeof(int) );
    gc_freeprop( snap2, "ID  " );

    /* shift all IDs of the second snapshot */
    for (i=snap1->nparttot; i<npartalltot; i++) {
      ids[i] += snap1->nparttot;
    }

    gadget_writeBlockInt( fout, "ID  ", npartalltot, 1, ids );

    free( ids );

    /* MASS and U */
    value = (float*)malloc( npartalltot * sizeof(float) );

    gc_readprop( snap1, "MASS", 0, snap1->nparttot );
    memcpy( value, gc_getdata( snap1 , "MASS" ), snap1->nparttot * sizeof(float) );
    gc_freeprop( snap1, "MASS" );

    gc_readprop( snap2, "MASS", 0, snap2->nparttot );
    memcpy( &value[snap1->nparttot], gc_getdata( snap2 , "MASS" ), snap2->nparttot * sizeof(float) );
    gc_freeprop( snap2, "MASS" );

    gadget_writeBlockFloat( fout, "MASS", npartalltot, 1, value );

    gc_readprop( snap1, "U   ", 0, snap1->nparttot );
    memcpy( value, gc_getdata( snap1 , "U   " ), snap1->nparttot * sizeof(float) );
    gc_freeprop( snap1, "U   " );

    gc_readprop( snap2, "U   ", 0, snap2->nparttot );
    memcpy( &value[snap1->nparttot], gc_getdata( snap2 , "U   " ), snap2->nparttot * sizeof(float) );
    gc_freeprop( snap2, "U   " );

    gadget_writeBlockFloat( fout, "U   ", npartalltot, 1, value );

    free( value );

    /* XNUC */
    value = (float*)malloc( nspecies * npartalltot * sizeof(float) );

    gc_readprop( snap1, "XNUC", 0, snap1->nparttot );
    memcpy( value, gc_getdata( snap1 , "XNUC" ), nspecies * snap1->nparttot * sizeof(float) );
    gc_freeprop( snap1, "XNUC" );

    gc_readprop( snap2, "XNUC", 0, snap2->nparttot );
    memcpy( &value[nspecies*snap1->nparttot], gc_getdata( snap2 , "XNUC" ), nspecies * snap2->nparttot * sizeof(float) );
    gc_freeprop( snap2, "XNUC" );

    gadget_writeBlockFloat( fout, "XNUC", npartalltot, nspecies, value );

    free( value );

    fclose( fout );
  } else {
    files = max( 1, snap1->multipleFiles ) + max( 1, snap2->multipleFiles );

    filenr = 0;

    filenr1 = 0;
    do {
      sprintf( fn, "%s.%d", fileout, filenr );
      fout = fopen( fn, "w" );

      gc_selectFile( snap1, filenr1 );
      printf( "Converting %s, snap %d into %s.\n", snap1->filename, filenr1, fn );

      gadget_writeHeader( fout, npartalltot, snap1->nparttot, files );

      gc_readprop( snap1, "POS ", 0, snap1->nparttot );
      value = (float*)gc_getdata( snap1, "POS " );
      if (off1) for (i=0; i<snap1->nparttot; i++) for (j=0; j<3; j++) value[i*3+j] = value[i*3+j]*scale1 + data_off1[j];
      gadget_writeBlockFloat( fout, "POS ", snap1->nparttot, 3, value );
      gc_freeprop( snap1, "POS " );

      gc_readprop( snap1, "VEL ", 0, snap1->nparttot );
      value = (float*)gc_getdata( snap1, "VEL " );
      if (voff1) for (i=0; i<snap1->nparttot; i++) for (j=0; j<3; j++) value[i*3+j] += data_voff1[j];
      gadget_writeBlockFloat( fout, "VEL ", snap1->nparttot, 3, value );
      gc_freeprop( snap1, "VEL " );

      gc_readprop( snap1, "ID  ", 0, snap1->nparttot );
      ids = (int*)gc_getdata( snap1, "ID  " );
      gadget_writeBlockInt( fout, "ID  ", snap1->nparttot, 1, ids );
      gc_freeprop( snap1, "ID  " );

      gc_readprop( snap1, "MASS", 0, snap1->nparttot );
      value = (float*)gc_getdata( snap1, "MASS" );
      gadget_writeBlockFloat( fout, "MASS", snap1->nparttot, 1, value );
      gc_freeprop( snap1, "MASS" );

      gc_readprop( snap1, "U   ", 0, snap1->nparttot );
      value = (float*)gc_getdata( snap1, "U   " );
      gadget_writeBlockFloat( fout, "U   ", snap1->nparttot, 1, value );
      gc_freeprop( snap1, "U   " );

      gc_readprop( snap1, "XNUC", 0, snap1->nparttot );
      value = (float*)gc_getdata( snap1, "XNUC" );
      gadget_writeBlockFloat( fout, "XNUC", snap1->nparttot, nspecies, value );
      gc_freeprop( snap1, "XNUC" );

      fclose( fout );

      filenr++;
      filenr1++;
    } while (filenr1 < snap1->multipleFiles);

    filenr2 = 0;
    do {
      sprintf( fn, "%s.%d", fileout, filenr );
      fout = fopen( fn, "w" );

      gc_selectFile( snap2, filenr2 );
      printf( "Converting %s, snap %d into %s.\n", snap2->filename, filenr2, fn );

      gadget_writeHeader( fout, npartalltot, snap2->nparttot, files );

      gc_readprop( snap2, "POS ", 0, snap2->nparttot );
      value = (float*)gc_getdata( snap2, "POS " );
      if (off2) for (i=0; i<snap2->nparttot; i++) for (j=0; j<3; j++) value[i*3+j] = value[i*3+j]*scale2 + data_off2[j];
      gadget_writeBlockFloat( fout, "POS ", snap2->nparttot, 3, value );
      gc_freeprop( snap2, "POS " );

      gc_readprop( snap2, "VEL ", 0, snap2->nparttot );
      value = (float*)gc_getdata( snap2, "VEL " );
      if (voff2) for (i=0; i<snap2->nparttot; i++) for (j=0; j<3; j++) value[i*3+j] += data_voff2[j];
      gadget_writeBlockFloat( fout, "VEL ", snap2->nparttot, 3, value );
      gc_freeprop( snap2, "VEL " );

      gc_readprop( snap2, "ID  ", 0, snap2->nparttot );
      ids = (int*)gc_getdata( snap2, "ID  " );
      for (i=0; i<snap2->nparttot; i++) ids[i] += snap1->npartalltot;
      gadget_writeBlockInt( fout, "ID  ", snap2->nparttot, 1, ids );
      gc_freeprop( snap2, "ID  " );

      gc_readprop( snap2, "MASS", 0, snap2->nparttot );
      value = (float*)gc_getdata( snap2, "MASS" );
      gadget_writeBlockFloat( fout, "MASS", snap2->nparttot, 1, value );
      gc_freeprop( snap2, "MASS" );

      gc_readprop( snap2, "U   ", 0, snap2->nparttot );
      value = (float*)gc_getdata( snap2, "U   " );
      gadget_writeBlockFloat( fout, "U   ", snap2->nparttot, 1, value );
      gc_freeprop( snap2, "U   " );

      gc_readprop( snap2, "XNUC", 0, snap2->nparttot );
      value = (float*)gc_getdata( snap2, "XNUC" );
      gadget_writeBlockFloat( fout, "XNUC", snap2->nparttot, nspecies, value );
      gc_freeprop( snap2, "XNUC" );

      fclose( fout );

      filenr++;
      filenr2++;
    } while (filenr2 < snap2->multipleFiles);
  }

  gc_deinit( snap1 );
  gc_deinit( snap2 );

  Py_RETURN_TRUE;
}

PyObject* _invertIDs(PyObject *self, PyObject *args) {
  PyArrayObject *pids;
  int *ids, *nids, *p, *pend, i;
  if (!PyArg_ParseTuple( args, "O!:invertIDs( ids )", &PyArray_Type, &pids )) {
    Py_RETURN_FALSE;
  }
  
  ids = (int*)PyArray_DATA(pids);

  nids = (int*)malloc( (PyArray_DIMS(pids)[0]+1) * sizeof(int) );
  nids[0] = -1;
  
  i = 0;
  pend = &ids[PyArray_DIMS(pids)[0]];
  for (p = ids; p != pend; p++) {
    nids[ *p ] = i;
    i++;
  }

  return (PyObject*)createPyIntArray( nids, PyArray_DIMS(pids)[0]+1 );
}

PyObject* _compareSnaps(PyObject *self, PyObject *args) {
  char *file1, *file2;
  t_gadgetSnap *snap1, *snap2;
  char pname[5];
  int prop, part, d, dim, match;
  float *data1, *data2;
  int *ids1, *ids2, *rids2, *p, *pend, i;

  if (!PyArg_ParseTuple( args, "ss:compareSnaps( filename1, filename2 )", &file1, &file2 )) {
    Py_RETURN_FALSE;
  }

  snap1 = gc_init( file1 );
  snap2 = gc_init( file2 );

  if (!snap1 || !snap2) Py_RETURN_FALSE;

  if (snap1->nparttot != snap2->nparttot) {
    printf( "Snapshot 1 and snapshot 2 contain different number of particles: %d | %d.\n", snap1->nparttot, snap2->nparttot );
    Py_RETURN_TRUE;
  }

  gc_readprop( snap1, "ID  ", 0, snap1->nparttot );
  ids1 = (int*)gc_getdata( snap1, "ID  " );
        gc_readprop( snap2, "ID  ", 0, snap2->nparttot );
  ids2 = (int*)gc_getdata( snap2, "ID  " );

  rids2 = (int*)malloc( (snap2->nparttot+1) * sizeof(int) );
  rids2[0] = -1;
  
  i = 0;
  pend = &ids2[snap2->nparttot];
  for (p = ids2; p != pend; p++) {
    rids2[ *p ] = i;
    i++;
  }

  for (prop=0; prop<gc_getpropcount(snap1); prop++) {
    gc_getpropname( snap1, prop, pname );
    if (!strcmp( pname, "HEAD" ) || !strcmp( pname, "ID  ")) {
      continue;
    }

    if (!gc_hasprop( snap2, pname )) {
      printf( "Property %s only present in snapshot 1.\n", pname );
      continue;
    }

    if (gc_getdim( snap1, pname ) != gc_getdim( snap2, pname )) {
      printf( "Property %s has different dimensions in both snapshots: %d | %d.\n", pname, gc_getdim( snap1, pname ), gc_getdim( snap2, pname ));
      continue;
    }

    gc_readprop( snap1, pname, 0, snap1->nparttot );
    data1 = (float*)gc_getdata( snap1, pname );
    gc_readprop( snap2, pname, 0, snap2->nparttot );
    data2 = (float*)gc_getdata( snap2, pname );
    
    match = 1;

    for (part=0; part<snap1->nparttot; part++) {
      dim = gc_getdim(snap1,pname);
      for (d=0; d<dim; d++) {
        if (data1[part*dim+d] != data2[rids2[ids1[part]]*dim+d]) {
          match = 0;
          break;
        }
      }
    }

    if (match) {
      printf( "Property %s is identical in both snapshots.\n", pname );
    } else {
      printf( "Property %s is different in both snapshots.\n", pname );
    }

    gc_freeprop( snap1, pname );
    gc_freeprop( snap2, pname );
  }

  for (prop=0; prop<gc_getpropcount(snap2); prop++) {
    gc_getpropname( snap2, prop, pname );

    if (!gc_hasprop( snap1, pname )) {
      printf( "Property %s only present in snapshot 2.\n", pname );
      continue;
    }
  }

  gc_deinit( snap1 );
  gc_deinit( snap2 );

  Py_RETURN_TRUE;
}

static PyMethodDef rgadgetMethods[] = {
  { "mergeICs", _mergeICs, METH_VARARGS, "" },
  { "invertIDs", _invertIDs, METH_VARARGS, "" },
  { "compareSnaps", _compareSnaps, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "rgadget", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  rgadgetMethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_rgadget(void)
{
  import_array();
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initrgadget(void)
{
  Py_InitModule( "rgadget", rgadgetMethods );
  import_array();
}
#endif
