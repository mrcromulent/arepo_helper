#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helm_eos.h"
#include "mersenne.h"
#include "pyhelm_eos.h"
#include "rgadget.h"
#include "pix2vec_ring.h"

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

typedef struct {
  double x, y, z, r;
} t_point;

int compare_points( const void* o1, const void* o2 ) {
  t_point *p1, *p2;
  p1 = (t_point*)o1;
  p2 = (t_point*)o2;

  if (p1->r < p2->r) return -1;
  else if (p1->r == p2->r) return 0;
  else return 1;
}

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

/* x[zone] is the left edge, x[zone-1] the right edge of a zone 
 *  => x[zone-1] <= px <= x[zone]
 * x[-1] is assumed to be 0
 * px is the real position */
float interpolate1D_valuedim( int zone, float px, int nzones, double *x, double *value, int valuedim, int valuenumber ) {
  float xcenter, xcenter2;
  int zone2;

  if (zone == 0) {
    xcenter = 0.5 * x[0];
  } else {
    xcenter = 0.5 * ( x[zone-1] + x[zone] );
  }

  if (px < xcenter) {
    if (zone == 0) return value[0];
    zone2 = zone-1;
  } else {
    if (zone == nzones-1) return value[nzones-1];
    zone2 = zone+1;
  }

  xcenter2 = 0.5 * ( x[zone2-1] + x[zone2] );
  
  return value[zone2*valuedim+valuenumber] + (px - xcenter2) / (xcenter - xcenter2) * (value[zone*valuedim+valuenumber] - value[zone2*valuedim+valuenumber]);
}

float interpolate1D( int zone, float px, int nzones, double *x, double *value ) {
  return interpolate1D_valuedim( zone, px, nzones, x, value, 1, 0 );
}

/* x and y contain the edges of the cells and are of length xsize+1 and ysize+2
 * if !reverse, array is assumed to be ordered as index = y*xsize+x, otherwise index = x*ysize+y */
float interpolate2D_valuedim( float px, float py, double *x, double *y, double *value, int xsize, int ysize, int reverse, int valuedim, int valuenumber ) {
  int i, ix1, ix2, iy1, iy2;
  int idx1, idx2, idx3, idx4;
  double val, val1, val2;
  double *cx, *cy;

  cx = (double*)malloc( xsize * sizeof(double) );
  cy = (double*)malloc( ysize * sizeof(double) );

  for (i=0; i<xsize; i++)
    cx[i] = 0.5 * ( x[i] + x[i+1] );
  for (i=0; i<ysize; i++)
    cy[i] = 0.5 * ( y[i] + y[i+1] );

  px = max( cx[0], min( cx[xsize-1], px ) );
  py = max( cy[0], min( cy[ysize-1], py ) );

  ix2 = 0;
  while (px > cx[ix2]) ix2++;
  ix1 = max( 0, ix2-1 );

  iy2 = 0;
  while (py > cy[iy2]) iy2++;
  iy1 = max( 0, iy2-1 );

  if (!reverse) {
    idx1 = (iy1*xsize + ix1)*valuedim + valuenumber;
    idx2 = (iy1*xsize + ix2)*valuedim + valuenumber;
    idx3 = (iy2*xsize + ix1)*valuedim + valuenumber;
    idx4 = (iy2*xsize + ix2)*valuedim + valuenumber;
  } else {
    idx1 = (ix1*ysize + iy1)*valuedim + valuenumber;
    idx2 = (ix1*ysize + iy2)*valuedim + valuenumber;
    idx3 = (ix2*ysize + iy1)*valuedim + valuenumber;
    idx4 = (ix2*ysize + iy2)*valuedim + valuenumber;
  }

  /* interpolate in x-direction first */
  if (ix1 == ix2) {
    val1 = value[idx1];
    val2 = value[idx3];
  } else {
    val1 = value[idx1] + (px - cx[ix1]) / (cx[ix2] - cx[ix1]) * (value[idx2] - value[idx1]);
    val2 = value[idx3] + (px - cx[ix1]) / (cx[ix2] - cx[ix1]) * (value[idx4] - value[idx3]);
  }

  /* interpolate in y-direction afterwards */
  if (iy1 == iy2) {
    val = val1;
  } else {
    val = val1 + (py - cy[iy1]) / (cy[iy2] - cy[iy1]) * (val2 - val1);
  }

  free( cx );
  free( cy );

  return val;
}

float interpolate2D( float px, float py, double *x, double *y, double *value, int xsize, int ysize, int reverse ) {
  return interpolate2D_valuedim( px, py, x, y, value, xsize, ysize, reverse, 1, 0 );
}

double* resize( double *data, int oldsize, int newsize ) {
  double *newdata;
  int copylength;
  
  newdata = (double*)malloc( newsize * sizeof(double) );
    
  copylength = min( oldsize, newsize );
  
  memcpy( newdata, data, copylength * sizeof(double) );
  if (newsize > oldsize) {
    memset( &newdata[oldsize], 0, (newsize-copylength)*sizeof(double) );
  }
    
  free(data);
  return newdata;
}

PyObject* _create_particles_cube(PyObject *self, PyObject *args) {
  PyObject *data, *dict;
  t_helm_eos_table *helm_eos_table;
  long long bs, i, ix, iy, iz;
  int npart, usecells, ncells, nspecies;
  double temp, pmass, mass, dmass;
  PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
  double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
  t_point *cube;
  int j, k, n, index;
  double bsh;
  double *ndata_pos, *ndata_mass, *ndata_u, *ndata_vel, *ndata_xnuc;
  int p, cellspercell;
  double last_radius, next_radius, lr3, nr3, radius, newradius;
  struct eos_result res;

  usecells = 0;
  nspecies = 3;
  temp = 0.0;
  pmass = 0.0;
  if (!PyArg_ParseTuple( args, "OO&i|iidd:create_particles_cube( data, eos, npart, [usecells, nspecies, temp, pmass] )", &data, &pyConvertHelmEos, &helm_eos_table, &npart, &usecells, &nspecies, &temp, &pmass )) {
    return 0;
  }

  data_u = data_H = data_HE = data_xnuc = NULL;
  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  rho = (PyArrayObject*)PyDict_GetItemString( data, "rho" );
  data_rho = (double*)PyArray_DATA(rho);

  if (temp == 0.0) {
    u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
    data_u = (double*)PyArray_DATA(u);
  }

  if (nspecies == 3) {
    H = (PyArrayObject*)PyDict_GetItemString( data, "H" );
    data_H = (double*)PyArray_DATA(H);
    HE = (PyArrayObject*)PyDict_GetItemString( data, "HE" );
    data_HE = (double*)PyArray_DATA(HE);
  } else {
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    data_xnuc = (double*)PyArray_DATA(xnuc);
  }

  mass = 0.0;
  for (i=0; i<ncells; i++) mass += data_dm[i];

  if (pmass > 0) {
    npart = floor( mass / pmass + 0.5 );
  } else {
    pmass = mass / npart;
  }

  bs = floor( pow( 3.0*npart, 1./3. ) );
  bsh = (bs - 1.0) / 2.0;
  cube = (t_point*)malloc( bs*bs*bs*sizeof( t_point ) );
  for (ix=0; ix<bs; ix++) {
    for (iy=0; iy<bs; iy++) {
      for (iz=0; iz<bs; iz++) {
        i = (iz*bs + iy)*bs + ix;
        cube[i].x = ix - bsh;
        cube[i].y = iy - bsh;
        cube[i].z = iz - bsh;
        cube[i].r = sqrt( cube[i].x*cube[i].x + cube[i].y*cube[i].y + cube[i].z*cube[i].z );
      }
    }
  }

  qsort( cube, bs*bs*bs, sizeof(t_point), compare_points );

  ndata_pos = (double*)malloc( 3*npart*2*sizeof(double) );
  ndata_mass = (double*)malloc( npart*2*sizeof(double) );
  ndata_u = (double*)malloc( npart*2*sizeof(double) );
  ndata_vel = (double*)malloc( 3*npart*2*sizeof(double) );
  ndata_xnuc = (double*)malloc( nspecies*npart*2*sizeof(double) );

  p = 0;
  last_radius = 0.0;
  
  if (!usecells || usecells > ncells) {
    usecells = ncells;
  }

  cellspercell = max( 1, floor( ncells / usecells + 0.5 ) );
  if (usecells && cellspercell > 1) {
    usecells = floor( ncells / cellspercell );
    printf( "Cells per cell: %d, using %d artificial cells.\n", cellspercell, usecells );
  }
  
  mass = 0; 
  for (i=0; i<usecells; i++) {
    for (j=i*cellspercell; j<(i+1)*cellspercell; j++)
      mass += data_dm[j];
    
    dmass = mass - p*pmass;
    n = floor( dmass / pmass + 0.5 );
  
    next_radius = data_r[(i+1)*cellspercell];
    lr3 = pow( last_radius, 3.0 );
    nr3 = pow( next_radius, 3.0 );
    
    index = i*cellspercell;
    for (j=0; j<n; j++) {
      radius = cube[p].r;
      newradius = pow( 1.0 * (j+1.) / n * (nr3 - lr3) + lr3, 1.0/3.0 );
      
      while ( data_r[index] < newradius )
        index++;
      
      if (radius > 0) {
        ndata_pos[p*3]   = cube[p].x / radius * newradius;
        ndata_pos[p*3+1] = cube[p].y / radius * newradius;
        ndata_pos[p*3+2] = cube[p].z / radius * newradius;
      } else {
        ndata_pos[p*3]   = cube[p].x;
        ndata_pos[p*3+1] = cube[p].y;
        ndata_pos[p*3+2] = cube[p].z;
      }
      
      ndata_mass[p] = pmass;
      if (nspecies == 3) {
        ndata_xnuc[p*3]   = data_H[index];
        ndata_xnuc[p*3+1] = data_HE[index];
        ndata_xnuc[p*3+2] = 1. - data_H[index] - data_HE[index];
      } else {
        for (k=0; k<nspecies; k++)
          ndata_xnuc[p*nspecies+k] = data_xnuc[k];
      }

      if (temp > 0.0) {
        eos_calc_tgiven( helm_eos_table, data_rho[index], &ndata_xnuc[p], temp, &res );
        ndata_u[p] = res.e.v;
      } else {
        ndata_u[p] = data_u[index];
      }
      
      p++;
    }
    last_radius = next_radius;
  }

  memset( ndata_vel, 0, 3 * p * sizeof(double) );

  free( cube );

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "pos", (PyObject*)createPyArray( ndata_pos, p, 3 ) );
  free( ndata_pos );
  PyDict_SetStolenItem( dict, "mass", (PyObject*)createPyArray( ndata_mass, p, 1 ) );
  free( ndata_mass );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( ndata_u, p, 1 ) );
  free( ndata_u );
  PyDict_SetStolenItem( dict, "vel", (PyObject*)createPyArray( ndata_vel, p, 3 ) );
  free( ndata_vel );
  PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( ndata_xnuc, p, nspecies ) );
  free( ndata_xnuc );
  PyDict_SetStolenItem( dict, "count", (PyObject*)PyLong_FromLong( p ) );
  
  return dict;
}

PyObject* _create_particles_healpix(PyObject *self, PyObject *args, PyObject *keywds) {
  PyObject *data, *dict;
  t_helm_eos_table *helm_eos_table;
  PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
  double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
  double pmass, mass, mtot, masslast, width;
  int npart, nspecies, ncells, np;
  double cx, cy, cz;
  double *ndata_pos, *ndata_mass, *ndata_rho, *ndata_u, *ndata_vel, *ndata_xnuc;
  int idxlast, i, j, k, p, found, index, npart_allocated;
  double n1, n2, n, nlast, rad1, rad2, rad, vec[3];
  int makebox, boxres, nbox, ix, iy, iz;
  double boxsize, hboxsize, boxcellsize, minenergy;
  int randomizeshells, randomizeradii, noxnuc, fixgridpressure, transition_done;
  double phi1, phi2, phi3, x, y, z, x2, y2, z2;
  double gridenergy, griddensity, gridtemp, radius, boxfactor;
  double radius_last, dr, vol, volboxcell, volcell;
  double transition_radius, transition_pmass;

  static char *kwlist[] = { "data", "eos", "npart", "nspecies", "pmass", "makebox", "boxsize", "boxres", "randomizeshells", "randomizeradii", "minenergy", "fixgridpressure", "griddensity", "boxfactor", "gridtemp", "cx", "cy", "cz", "transition_radius", "transition_pmass", NULL };
  
  nspecies = 3;
  pmass = 0;
  makebox = 0;
  boxsize = 0;
  boxres = 16;
  randomizeshells = 0;
  randomizeradii = 0;
  minenergy = 0;
  fixgridpressure = 0;
  griddensity = 1e-5;
  boxfactor = 10;
  gridtemp = 1e8;
  cx = cy = cz = 0;
  transition_radius = 0;
  transition_pmass = 0;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "OO&i|ididiiididddddddd:create_particles_healpix( data, eos, npart, [nspecies, pmass, makebox, boxsize, boxres, randomizeshells, randomizeradii, minenergy, fixgridpressure, griddensity, boxfactor, gridtemp, cx, cy, cz, transition_radius, transition_pmass] )", kwlist, &data, &pyConvertHelmEos, &helm_eos_table, &npart, &nspecies, &pmass, &makebox, &boxsize, &boxres, &randomizeshells, &randomizeradii, &minenergy, &fixgridpressure, &griddensity, &boxfactor, &gridtemp, &cx, &cy, &cz, &transition_radius, &transition_pmass )) {
    return 0;
  }
  
  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  rho = (PyArrayObject*)PyDict_GetItemString( data, "rho" );
  data_rho = (double*)PyArray_DATA(rho);
  u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
  data_u = (double*)PyArray_DATA(u);

  data_H = data_HE = data_xnuc = NULL;
  if (nspecies == 3 && !PyDict_Contains( data, PyUnicode_FromString( "xnuc" ))) {
    noxnuc = 1;
    H = (PyArrayObject*)PyDict_GetItemString( data, "H" );
    data_H = (double*)PyArray_DATA(H);
    HE = (PyArrayObject*)PyDict_GetItemString( data, "HE" );
    data_HE = (double*)PyArray_DATA(HE);
  } else {
    noxnuc = 0;
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    data_xnuc = (double*)PyArray_DATA(xnuc);
  }

  mtot = 0.0;
  for (i=0; i<ncells; i++) mtot += data_dm[i];

  printf( "Total mass: %g solar masses\n", mtot / 1.989e33 );

  if (pmass > 0) {
    npart = floor( mtot / pmass + 0.5 );
  } else {
    pmass = mtot / npart;
  }
  
  if (transition_radius > 0) {
          transition_done = 0;
    printf( "Doing transition at r=%g to pmass=%g\n", transition_radius, transition_pmass );
  } else {
    transition_done = 1;
    printf( "Using %d particles of mass %g (%g solar masses)\n", npart, pmass, pmass/1.989e33 );
  }

  seed();
  
  npart_allocated = npart;
  ndata_pos = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_mass = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_rho = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_u = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_vel = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_xnuc = (double*)malloc( nspecies*npart_allocated*sizeof(double) );
  p = 0;
  
  idxlast = 0;
  nlast = 0;
  masslast = 0;
  mass = 0;
  rad1 = data_r[0];
  rad2 = rad1;
  found = 0;
  volcell = 0;
  phi1 = phi2 = phi3 = 0;
  np = 0;
  
  i = 0;
  while (i < ncells-1) {
    rad2 = data_r[i+1]; /* increases with index */
    n1 = sqrt(M_PI/12.) * (rad2+rad1) / (rad2-rad1); /* decreases with index */
    
    if ((!transition_done) && (rad2 > transition_radius)) {
      pmass = transition_pmass;
      transition_done = 1;
      printf( "Transition done at rad2=%g, index=%d, mass=%g, particles before: %d.\n", rad2, i, mass/1.989e33, p);
    }
    
    mass += data_dm[i]; /* increases with index */
    n2 = sqrt(mass/pmass/12.); /* increases with index */
    
    if ( floor( n2 ) > nlast ) {
      nlast = floor( n2 );
      idxlast = i;
      masslast = mass;
    }
    
    if ( n2 > n1 && !found ) {
      n = floor( n2 + 0.5 );
      found = 1;
    }
    
    if (found && nlast >= n) {
      i = idxlast;
      mass = masslast;
      rad2 = data_r[i+1];
      
      rad = 0.5 * (rad2+rad1);
      width = rad2 - rad1;
      index = i;
      while (data_r[index] > rad) index--;
      
      np = 12 * n * n;
      
      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }
      
      while (p + np >= npart_allocated) {
        npart_allocated *= 2;

        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_rho = (double*)realloc( ndata_rho, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }
      
      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );
        
        double prad = rad;
        if (randomizeradii)
                  prad += 0.1 * width * (randMT() - 0.5);

        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
          
          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
          ndata_pos[p*3+0] = prad*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[p*3+1] = prad*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[p*3+2] = prad*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[p*3+k] = vec[k] * prad;
        }
        
        if(makebox) {
          int inbox = 1;
          for (k=0; k<3; k++)
            if(ndata_pos[p*3+k] < -0.5*boxsize || ndata_pos[p*3+k] > +0.5*boxsize) {
              inbox = 0;
              break;
             }
          if(!inbox)
            continue;
        }
        
        ndata_mass[p] = pmass;
        ndata_rho[p] = data_rho[index];
        ndata_u[p] = data_u[index];
        if (noxnuc) {
          ndata_xnuc[p*3]   = data_H[index];
          ndata_xnuc[p*3+1] = data_HE[index];
          ndata_xnuc[p*3+2] = 1. - data_H[index] - data_HE[index];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[p*nspecies+k] = data_xnuc[k];
        }
        
        p++;
      }
      
      volcell = 4./3. * M_PI * (rad2*rad2*rad2 - rad1*rad1*rad1) / np;
    
      rad1 = rad2;
      mass -= pmass * 12. * n * n;
      found = 0;
      nlast = 0;
    }
    
    i++;
  }
  
  n = floor( sqrt(mass/pmass/12.) + 0.5 );
  /* put last point near radius of star 
  if (ncells > 50) {
    rad1 = data_r[i - 6];
  } else {
    rad1 = data_r[i - 1];
  }
  */
  rad = 0.5 * (rad2+rad1);
  width = rad2 - rad1;
  index = i;
  while (data_r[index] > rad) index--;
  
  if (randomizeshells) {
    phi1 = randMT() * 2. * M_PI;
    phi2 = randMT() * 2. * M_PI;
    phi3 = randMT() * 2. * M_PI;
  }
  
  double prad = rad;
  if (randomizeradii)
                prad += 0.1 * width * (randMT() - 0.5);

  while (p + np >= npart_allocated) {
    npart_allocated *= 2;

    ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
    ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
    ndata_rho = (double*)realloc( ndata_rho, npart_allocated*sizeof(double) );
    ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
    ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
    ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
  }
  
  np = 12 * n * n;
  for (j=0; j<np; j++) {
    pix2vec_ring( n, j, vec );
    
    if (randomizeshells) {
      x = vec[0];
      y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
      z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
          
      x2 =  cos( phi2 )*x + sin( phi2 )*z;
      y2 = y;
      z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
      ndata_pos[p*3+0] = prad*( cos( phi3 )*x2 - sin( phi3 )*y2 );
      ndata_pos[p*3+1] = prad*( sin( phi3 )*x2 + cos( phi3 )*y2 );
      ndata_pos[p*3+2] = prad*z2;
    } else {    
      for (k=0; k<3; k++)
        ndata_pos[p*3+k] = vec[k] * prad;
    }
    
    if(makebox) {
      int inbox = 1;
      for (k=0; k<3; k++)
      if(ndata_pos[p*3+k] < -0.5*boxsize || ndata_pos[p*3+k] > +0.5*boxsize) {
        inbox = 0;
        break;
      }
      if(inbox == 0)
        continue;
    }
    
    ndata_mass[p] = pmass;
    ndata_rho[p] = data_rho[index];
    
    ndata_u[p] = data_u[index];
    if (ndata_u[p] < minenergy) ndata_u[p] = minenergy;
    if (noxnuc) {
      ndata_xnuc[p*3]   = data_H[index];
      ndata_xnuc[p*3+1] = data_HE[index];
      ndata_xnuc[p*3+2] = 1. - data_H[index] - data_HE[index];
    } else {
      for (k=0; k<nspecies; k++)
        ndata_xnuc[p*nspecies+k] = data_xnuc[k];
    }
    
    p++;
  }
  
  if (makebox) {
    double *grid_xnuc = malloc( nspecies * sizeof(double) );
    if (noxnuc) {
            grid_xnuc[0] = data_H[index];
      grid_xnuc[1] = data_HE[index];
      grid_xnuc[2] = 1. - data_H[index] - data_HE[index];
    } else {
            for (k=0; k<nspecies; k++)
              grid_xnuc[k] = data_xnuc[k];
    }
    
                if (fixgridpressure) {
            struct eos_result res;
      double pressure, temp;

      temp = -1.;
      eos_calc_egiven( helm_eos_table, data_rho[index], grid_xnuc, data_u[index], &temp, &res );
      printf( "rho=%g, u=%g, temp=%g, p=%g\n", data_rho[index], data_u[index], temp, res.p.v );

      temp = -1;
      pressure = res.p.v;
      eos_calc_pgiven( helm_eos_table, griddensity, grid_xnuc, pressure, &temp, &res);
      printf( "rho=%g, u=%g, temp=%g, p=%g\n", griddensity, res.e.v, temp, res.p.v );
      gridenergy = res.e.v;

      if (temp == -1)
                                temp = gridtemp;
      eos_calc_tgiven( helm_eos_table, griddensity, grid_xnuc, temp, &res );
      printf( "rho=%g, u=%g, temp=%g, p=%g\n", griddensity, res.e.v, temp, res.p.v );
      gridenergy = res.e.v;
    } else {
      gridenergy = minenergy;
    }
    
    hboxsize = boxsize / 2.;
    boxcellsize = boxsize / (double)boxres;
    volboxcell = boxcellsize * boxcellsize * boxcellsize;

    radius_last = rad2;
    //volcell = 4./3. * M_PI * (rad2*rad2*rad2 - rad1*rad1*rad1);
    printf( "volcell=%g, volboxcell=%g\n", volcell, volboxcell );
    volcell *= boxfactor;
    while (volcell < volboxcell && boxfactor > 0.)
    {
      dr = pow( volcell, 1./3. );

      vol = 4./3. * M_PI * (pow(radius_last+dr,3)-pow(radius_last,3));
      n = floor( sqrt( vol / (12. * volcell) ) + 0.5 );
      np = 12 * n * n;
      radius = radius_last + 0.5 * dr;

      printf( "volcell=%g, volboxcell=%g, np=%d\n", volcell, volboxcell, np );

      while (p + np >= npart_allocated) {
        npart_allocated *= 2;

        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_rho = (double*)realloc( ndata_rho, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }

      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }

      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );

        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];

          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;

          ndata_pos[p*3+0] = radius*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[p*3+1] = radius*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[p*3+2] = radius*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[p*3+k] = vec[k] * radius;
        }

                                int inbox = 1;
                                for (k=0; k<3; k++)
                                        if(ndata_pos[p*3+k] < -hboxsize || ndata_pos[p*3+k] > +hboxsize) {
                                                inbox = 0;
                                                break;
                                         }
                                if(inbox == 0)
                                        continue;
                                
        ndata_rho[p] = griddensity;
        ndata_mass[p] = griddensity * volcell;
        
        for (k=0; k<nspecies; k++)
          ndata_xnuc[p*nspecies+k] = grid_xnuc[k];

        p++;
      }

      radius_last += dr;
      volcell *= boxfactor;
    }
    
    for (i=0; i<p; i++) {
      ndata_pos[i*3+0] += hboxsize + cx;
      ndata_pos[i*3+1] += hboxsize + cy;
      ndata_pos[i*3+2] += hboxsize + cz;
    }

    nbox = boxres * boxres * boxres;

    while (p + nbox >= npart_allocated) {
      npart_allocated *= 2;

      ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
      ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
      ndata_rho = (double*)realloc( ndata_rho, npart_allocated*sizeof(double) );
      ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
      ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
      ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
    }

    char *box = malloc( nbox );
    memset( box, 0, nbox );
    for (i=0; i<p; i++) {
      ix = floor( ndata_pos[i*3+0] / boxcellsize );
      iy = floor( ndata_pos[i*3+1] / boxcellsize );
      iz = floor( ndata_pos[i*3+2] / boxcellsize );

      box[iz*boxres*boxres + iy*boxres + ix] = 1;
    }

    int boxCount = 0;
    for (i=0; i<nbox; i++) {
      ix = i % boxres;
      iy = (i / boxres) % boxres;
      iz = i / (boxres*boxres);
      
      x = (ix+0.5) * boxcellsize;
      y = (iy+0.5) * boxcellsize;
      z = (iz+0.5) * boxcellsize;
      
      if(!box[i]) {
        ndata_pos[p*3+0] = x;
        ndata_pos[p*3+1] = y;
        ndata_pos[p*3+2] = z;
        
        ndata_rho[p] = griddensity;
        ndata_mass[p] = griddensity * volboxcell;
        
        ndata_u[p] = gridenergy;
        for (k=0; k<nspecies; k++)
          ndata_xnuc[p*nspecies+k] = grid_xnuc[k];
        
        p++;
        boxCount++;
      }
    }
    free(box);
    printf( "Added %d box cells.\n", boxCount );
    
    free( grid_xnuc );
  }
  
  memset( ndata_vel, 0, 3 * p * sizeof(double) );
  
  printf( "Created %d particles.\n", p );

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "pos", (PyObject*)createPyArray( ndata_pos, p, 3 ) );
  free( ndata_pos );
  PyDict_SetStolenItem( dict, "mass", (PyObject*)createPyArray( ndata_mass, p, 1 ) );
  free( ndata_mass );
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( ndata_rho, p, 1 ) );
  free( ndata_rho );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( ndata_u, p, 1 ) );
  free( ndata_u );
  PyDict_SetStolenItem( dict, "vel", (PyObject*)createPyArray( ndata_vel, p, 3 ) );
  free( ndata_vel );
  PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( ndata_xnuc, p, nspecies ) );
  free( ndata_xnuc );
  PyDict_SetStolenItem( dict, "count", (PyObject*)PyLong_FromLong( p ) );
  
  return dict;
}

PyObject* _create_particles_healpix_vol(PyObject *self, PyObject *args, PyObject *keywds) {
  PyObject *data, *dict;
  t_helm_eos_table *helm_eos_table;
  PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
  double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
  double *ndata_pos, *ndata_mass, *ndata_u, *ndata_vel, *ndata_xnuc;
  int ncells, nspecies, npart_allocated, np, npart_guess;
  double maxrad, vol, voltot, volshell, cellwidth, n;
  double mtot, cx, cy, cz;
  double rad_inner, rad_outer, rad, dr, vec[3];
  int makebox, boxres, nbox, ix, iy, iz, i, j, k, p, index;
  double boxsize, hboxsize, boxcellsize, minenergy;
  int randomizeshells, randomizeradii, noxnuc;
  double phi1, phi2, phi3, x, y, z, x2, y2, z2;
  double gridenergy, griddensity, gridtemp, radius, boxfactor;
  double radius_last, volboxcell, volcell;

  static char *kwlist[] = { "data", "eos", "vol", "nspecies", "makebox", "boxsize", "boxres", "randomizeshells", "randomizeradii", "minenergy", "griddensity", "boxfactor", "gridtemp", "cx", "cy", "cz", NULL };
  
  nspecies = 3;
  makebox = 0;
  boxsize = 0;
  boxres = 16;
  randomizeshells = 0;
  randomizeradii = 0;
  minenergy = 0;
  griddensity = 1e-5;
  boxfactor = 10;
  gridtemp = 1e8;
  cx = cy = cz = 0;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "OO&d|iidiiiddddddd:create_particles_healpix_vol( data, eos, vol, [nspecies, makebox, boxsize, boxres, randomizeshells, randomizeradii, minenergy, griddensity, boxfactor, gridtemp, cx, cy, cz] )", kwlist, &data, &pyConvertHelmEos, &helm_eos_table, &vol, &nspecies, &makebox, &boxsize, &boxres, &randomizeshells, &randomizeradii, &minenergy, &griddensity, &boxfactor, &gridtemp, &cx, &cy, &cz )) {
    return 0;
  }
  
  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  rho = (PyArrayObject*)PyDict_GetItemString( data, "rho" );
  data_rho = (double*)PyArray_DATA(rho);
  u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
  data_u = (double*)PyArray_DATA(u);

  data_H = data_HE = data_xnuc = NULL;
  if (nspecies == 3 && !PyDict_Contains( data, PyUnicode_FromString( "xnuc" ))) {
    noxnuc = 1;
    H = (PyArrayObject*)PyDict_GetItemString( data, "H" );
    data_H = (double*)PyArray_DATA(H);
    HE = (PyArrayObject*)PyDict_GetItemString( data, "HE" );
    data_HE = (double*)PyArray_DATA(HE);
  } else {
    noxnuc = 0;
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    data_xnuc = (double*)PyArray_DATA(xnuc);
  }

  mtot = 0.0;
  for (i=0; i<ncells; i++) mtot += data_dm[i];

  printf( "Total mass: %g solar masses\n", mtot / 1.989e33 );
  
  seed();

  maxrad = data_r[ncells-1];
  voltot = 4. / 3. * M_PI * maxrad * maxrad * maxrad;
  npart_guess = (int)floor(voltot/vol+0.5);

  printf( "Volume per cell: %g, Total Volume: %g, Guess for number of cells: %d\n", vol, voltot, npart_guess );
  
  cellwidth = pow( vol, 1./3. );
  
  npart_allocated = npart_guess;
  ndata_pos = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_mass = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_u = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_vel = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_xnuc = (double*)malloc( nspecies*npart_allocated*sizeof(double) );
  p = 0;
  
  rad_inner = 0;
  rad_outer = 0;
  
  index = 0;
  phi1 = phi2 = phi3 = 0;
  np = 0;
  
  i = 0;
  while (i < ncells-1) {
    rad_outer = data_r[i+1];
    dr = rad_outer - rad_inner;
    
    if (dr > cellwidth || i == ncells-1) {
      rad = 0.5 * (rad_inner+rad_outer);
      index = i;
      while (data_r[index] < rad) index--;
      
      volshell = 4. / 3. * M_PI * ( rad_outer*rad_outer*rad_outer - rad_inner*rad_inner*rad_inner );
      n = floor( sqrt( volshell / (12. * vol) ) + 0.5 );
      np = 12 * n * n;

      printf( "Adding shell at r=%g with %d particles.\n", rad, np );
      
      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }
      
      while (p + np >= npart_allocated) {
        npart_allocated *= 2;

        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }
      
      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );
        
        double prad = rad;
        if (randomizeradii) {
                prad += 0.2 * cellwidth * (randMT() - 0.5);
          index = i;
          while (data_r[index] < rad) index--;
        }

        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
          
          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
          ndata_pos[p*3+0] = prad*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[p*3+1] = prad*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[p*3+2] = prad*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[p*3+k] = vec[k] * prad;
        }
        
        ndata_mass[p] = data_rho[index];
        ndata_u[p] = data_u[index];
        if (noxnuc) {
          ndata_xnuc[p*3]   = data_H[index];
          ndata_xnuc[p*3+1] = data_HE[index];
          ndata_xnuc[p*3+2] = 1. - data_H[index] - data_HE[index];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[p*nspecies+k] = data_xnuc[k];
        }
        
        p++;
      }
    
      rad_inner = rad_outer;
    }
    
    i++;
  }

  if (makebox) {
    double *grid_xnuc = malloc( nspecies * sizeof(double) );
    if (noxnuc) {
      grid_xnuc[0] = data_H[index];
      grid_xnuc[1] = data_HE[index];
      grid_xnuc[2] = 1. - data_H[index] - data_HE[index];
    } else {
      for (k=0; k<nspecies; k++)
        grid_xnuc[k] = data_xnuc[k];
    }
          
    gridenergy = minenergy;
    
    hboxsize = boxsize / 2.;
    boxcellsize = boxsize / (double)boxres;
    volboxcell = boxcellsize * boxcellsize * boxcellsize;

    radius_last = rad_outer;
    volcell = vol;
    printf( "volcell=%g, volboxcell=%g\n", volcell, volboxcell );
    volcell *= boxfactor;
    while (volcell < volboxcell)
    {
      dr = pow( volcell, 1./3. );

      vol = 4./3. * M_PI * (pow(radius_last+dr,3)-pow(radius_last,3));
      n = floor( sqrt( vol / (12. * volcell) ) + 0.5 );
      np = 12 * n * n;
      radius = radius_last + 0.5 * dr;

      printf( "volcell=%g, volboxcell=%g, np=%d\n", volcell, volboxcell, np );

      while (p + np >= npart_allocated) {
        npart_allocated *= 2;

        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }

      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }

      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );

        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];

          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;

          ndata_pos[p*3+0] = radius*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[p*3+1] = radius*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[p*3+2] = radius*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[p*3+k] = vec[k] * radius;
        }

        ndata_mass[p] = griddensity;
        
        for (k=0; k<nspecies; k++)
          ndata_xnuc[p*nspecies+k] = grid_xnuc[k];

        p++;
      }

      radius_last += dr;
      volcell *= boxfactor;
    }
    
    for (i=0; i<p; i++) {
      ndata_pos[i*3+0] += hboxsize + cx;
      ndata_pos[i*3+1] += hboxsize + cy;
      ndata_pos[i*3+2] += hboxsize + cz;
    }

    nbox = boxres * boxres * boxres;
    
    for (i=0; i<nbox; i++) {
            ix = i % boxres;
      iy = (i / boxres) % boxres;
      iz = i / (boxres*boxres);
        
      x = (ix+0.5) * boxcellsize - hboxsize;
      y = (iy+0.5) * boxcellsize - hboxsize;
      z = (iz+0.5) * boxcellsize - hboxsize;
      
      radius = sqrt( x*x + y*y + z*z );

      if (radius > radius_last) {
        while (p + np >= npart_allocated) {
          npart_allocated *= 2;

          ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
          ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
          ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
          ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
          ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
        }

        ndata_pos[p*3+0] = x + hboxsize;
        ndata_pos[p*3+1] = y + hboxsize;
        ndata_pos[p*3+2] = z + hboxsize;
        
        ndata_mass[p] = griddensity;
        
        ndata_u[p] = gridenergy;
        for (k=0; k<nspecies; k++)
          ndata_xnuc[p*nspecies+k] = grid_xnuc[k];
        
        p++;
        
      }
    }
    
    free( grid_xnuc );
  }
  
  memset( ndata_vel, 0, 3 * p * sizeof(double) );
  
  printf( "Created %d particles.\n", p );

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "pos", (PyObject*)createPyArray( ndata_pos, p, 3 ) );
  free( ndata_pos );
  PyDict_SetStolenItem( dict, "mass", (PyObject*)createPyArray( ndata_mass, p, 1 ) );
  free( ndata_mass );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( ndata_u, p, 1 ) );
  free( ndata_u );
  PyDict_SetStolenItem( dict, "vel", (PyObject*)createPyArray( ndata_vel, p, 3 ) );
  free( ndata_vel );
  PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( ndata_xnuc, p, nspecies ) );
  free( ndata_xnuc );
  PyDict_SetStolenItem( dict, "count", (PyObject*)PyLong_FromLong( p ) );
  
  return dict;
}

PyObject* _create_particles_healpix_grad(PyObject *self, PyObject *args, PyObject *keywds) {
  PyObject *data, *dict;
  t_helm_eos_table *helm_eos_table;
  PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
  double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
  double cellsperdecade, gradlimit;
  int nspecies, makebox, boxres, randomizeshells;
  double boxsize, boxfactor, minenergy, minvolume, maxmass, roffset;
  int ncells, cell, noxnuc, nshells;
  int npart_allocated, npart, *box;
  double *ndata_pos, *ndata_mass, *ndata_u, *ndata_vel, *ndata_xnuc;
  double radius_last, dr, rholast, ulast, drho, du, dmass;
  int n, np, i, j, k, index;
  double vol, volcell, radius, phi1, phi2, phi3;
  double x, y, z, x2, y2, z2, vec[3];
  int nbox, idx, ix, iy, iz;
  double px, py, pz, cx, cy, cz;
  double hboxsize, boxcellsize, mtot, volboxcell, vollast, width;

  static char *kwlist[] = { "data", "eos", "nspecies", "cellsperdecade", "makebox", "boxsize", "boxres", "boxfactor", "randomizeshells", "minenergy", "minvolume", "maxmass", "roffset", NULL };
  
  nspecies = 3;
  cellsperdecade = 10;
  makebox = 0;
  boxsize = 0;
  boxres = 16;
  boxfactor = 10.;
  randomizeshells = 0;
  minenergy = 0;
  minvolume = 0;
  maxmass = 0;
  roffset = 0;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "OO&|ididididddd:create_particles_healpix_grad( data, eos, [nspecies, cellsperdecade, makebox, boxsize, boxres, boxfactor, randomizeshells, minenergy, minvolume, maxmass, roffset] )", kwlist, &data, &pyConvertHelmEos, &helm_eos_table, &nspecies, &cellsperdecade, &makebox, &boxsize, &boxres, &boxfactor, &randomizeshells, &minenergy, &minvolume, &maxmass, &roffset )) {
    return 0;
  }
  
  gradlimit = 1. / cellsperdecade;
  
  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  rho = (PyArrayObject*)PyDict_GetItemString( data, "rho" );
  data_rho = (double*)PyArray_DATA(rho);
  u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
  data_u = (double*)PyArray_DATA(u);

  data_H = data_HE = data_xnuc = NULL;
  if (nspecies == 3 && !PyDict_Contains( data, PyUnicode_FromString( "xnuc" ))) {
    noxnuc = 1;
    H = (PyArrayObject*)PyDict_GetItemString( data, "H" );
    data_H = (double*)PyArray_DATA(H);
    HE = (PyArrayObject*)PyDict_GetItemString( data, "HE" );
    data_HE = (double*)PyArray_DATA(HE);
  } else {
    noxnuc = 0;
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    data_xnuc = (double*)PyArray_DATA(xnuc);
  }

  mtot = 0.0;
  for (i=0; i<ncells; i++) mtot += data_dm[i];

  printf( "Total mass: %g solar masses\n", mtot / 1.989e33 );
  printf( "Radius: %g cm\n", data_r[ncells-1] );

  npart_allocated = 10000;
  ndata_pos = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_mass = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_u = (double*)malloc( npart_allocated*sizeof(double) );
  ndata_vel = (double*)malloc( 3*npart_allocated*sizeof(double) );
  ndata_xnuc = (double*)malloc( nspecies*npart_allocated*sizeof(double) );
  
  npart = 0;
  radius = 0;
  radius_last = 0;
  nshells = 0;
  phi1 = phi2 = phi3 = 0;
  
  dmass = data_dm[0];
  rholast = data_rho[0];
  ulast = data_u[0];
    
  cell = 0;
  vollast = 0;
  while (cell < ncells-1) {
    cell++;
    drho = (rholast - data_rho[cell]) / rholast;
    du = (ulast - data_u[cell]) / ulast;
    dmass += data_dm[cell];
    
    int do_shell = 0;

    if (fabs(drho) > gradlimit || fabs(du) > gradlimit) 
      do_shell = 1;

    if (maxmass > 0) {
      double size = data_r[cell] - radius_last;
      vol = 4./3. * M_PI * (pow(data_r[cell],3)-pow(radius_last,3));
      volcell = size*size*size;
      double cellmass = dmass * volcell / vol;

      if (cellmass > 0.9 * maxmass)
        do_shell = 2;
    }

    if (do_shell) {
      /* create a shell */
      n = floor( 0.5 + sqrt(M_PI/3.) * radius / (data_r[cell] - radius_last));
      np = 12 * n * n;

      radius = 0.5 * (data_r[cell] + radius_last);
      width = data_r[cell] - radius_last;
      vol = 4./3. * M_PI * (pow(data_r[cell],3)-pow(radius_last,3));
      
      volcell = vol / np;
      if (volcell < minvolume)
        continue;
      
      while (npart + np >= npart_allocated) {
        npart_allocated *= 2;
        
        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }
      
      index = cell;
      while (data_r[index] > radius) index--;

      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }

      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );

        double rad = radius;
        if (roffset > 0) {
          rad = radius + (randMT() - 0.5) * roffset * width;

          index = cell;
          while (data_r[index] > rad) index--;
        }
        
        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
          
          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
          ndata_pos[npart*3+0] = rad*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[npart*3+1] = rad*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[npart*3+2] = rad*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[npart*3+k] = vec[k] * radius;
        }
        
        ndata_mass[npart] = data_rho[index];
        
        ndata_u[npart] = max( data_u[index], minenergy );
        if (noxnuc) {
          ndata_xnuc[npart*3]   = data_H[index];
          ndata_xnuc[npart*3+1] = data_HE[index];
          ndata_xnuc[npart*3+2] = 1. - data_H[index] - data_HE[index];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[npart*nspecies+k] = data_xnuc[k];
        }
        
        npart++;
      }
      
      nshells++;
      printf( "shell %d (%d): radius=%g, dr=%g, np=%d, density=%g, u=%g, volume per cell=%g/%g\n", nshells, do_shell, radius, data_r[cell] - radius_last, np, data_rho[index], data_u[index], volcell, pow(data_r[cell]-radius_last,3) );
      vollast = volcell;
      
      dmass = 0;
      radius_last = data_r[cell];
      rholast = data_rho[cell];
      ulast = data_u[cell];
    }
  }
  
  volcell = vollast;
  if (makebox) {
    hboxsize = boxsize / 2.;
    boxcellsize = boxsize / (double)boxres;
    volboxcell = boxcellsize * boxcellsize * boxcellsize;
    
    volcell *= boxfactor;
    while (volcell < volboxcell)
    {
      dr = pow( volcell, 1./3. );
      
      vol = 4./3. * M_PI * (pow(radius_last+dr,3)-pow(radius_last,3));
      n = floor( sqrt( vol / (12. * volcell) ) + 0.5 );
      np = 12 * n * n;
      radius = radius_last + 0.5 * dr;
      
      printf( "volcell=%g, volboxcell=%g, np=%d\n", volcell, volboxcell, np );
      
      while (npart + np >= npart_allocated) {
        npart_allocated *= 2;
        
        ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
        ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
        ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
        ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
        ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
      }
      
      if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
      }
      
      for (j=0; j<np; j++) {
        pix2vec_ring( n, j, vec );
        
        if (randomizeshells) {
          x = vec[0];
          y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
          z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
          
          x2 =  cos( phi2 )*x + sin( phi2 )*z;
          y2 = y;
          z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
          ndata_pos[npart*3+0] = radius*( cos( phi3 )*x2 - sin( phi3 )*y2 );
          ndata_pos[npart*3+1] = radius*( sin( phi3 )*x2 + cos( phi3 )*y2 );
          ndata_pos[npart*3+2] = radius*z2;
        } else {
          for (k=0; k<3; k++)
            ndata_pos[npart*3+k] = vec[k] * radius;
        }
        
        ndata_mass[npart] = 0.;
        
        ndata_u[npart] = minenergy;
        if (noxnuc) {
          ndata_xnuc[npart*3]   = data_H[ncells-1];
          ndata_xnuc[npart*3+1] = data_HE[ncells-1];
          ndata_xnuc[npart*3+2] = 1. - data_H[ncells-1] - data_HE[ncells-1];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[npart*nspecies+k] = data_xnuc[k];
        }
        
        npart++;
      }
      
      radius_last += dr;
      volcell *= boxfactor;
    }
    
    nbox = boxres * boxres * boxres;
    
    box = (int*)malloc( nbox * sizeof(int) );
    memset( box, 0, nbox * sizeof(int) );
    
    for (i=0; i<npart; i++) {
      ndata_pos[i*3+0] += hboxsize;
      ndata_pos[i*3+1] += hboxsize;
      ndata_pos[i*3+2] += hboxsize;
      
      px = ndata_pos[i*3+0];
      py = ndata_pos[i*3+1];
      pz = ndata_pos[i*3+2];
      
      ix = min( max( floor( px / boxcellsize ), 0 ), boxres-1 );
      iy = min( max( floor( py / boxcellsize ), 0 ), boxres-1 );
      iz = min( max( floor( pz / boxcellsize ), 0 ), boxres-1 );
      
      cx = (ix+0.5) * boxcellsize;
      cy = (iy+0.5) * boxcellsize;
      cz = (iz+0.5) * boxcellsize;

      idx = 0;
      if (px > cx) idx += 1;
      if (py > cy) idx += 2;
      if (pz > cz) idx += 4;

      box[(boxres*iz+iy)*boxres+ix] |= 1 << idx;
    }
    
    for (i=0; i<nbox; i++) {
      if (box[i] != 0xFF) {
        ix = i % boxres;
        iy = (i / boxres) % boxres;
        iz = i / (boxres*boxres);
        
        if (npart + 1 >= npart_allocated) {
          npart_allocated *= 2;

          ndata_pos = (double*)realloc( ndata_pos, 3*npart_allocated*sizeof(double) );
          ndata_mass = (double*)realloc( ndata_mass, npart_allocated*sizeof(double) );
          ndata_u = (double*)realloc( ndata_u, npart_allocated*sizeof(double) );
          ndata_vel = (double*)realloc( ndata_vel, 3*npart_allocated*sizeof(double) );
          ndata_xnuc = (double*)realloc( ndata_xnuc, nspecies*npart_allocated*sizeof(double) );
        }
        
        ndata_pos[npart*3+0] = (ix+0.5) * boxcellsize;
        ndata_pos[npart*3+1] = (iy+0.5) * boxcellsize;
        ndata_pos[npart*3+2] = (iz+0.5) * boxcellsize;
        
        ndata_mass[npart] = 0;
        ndata_u[npart] = minenergy;
        if (noxnuc) {
          ndata_xnuc[npart*3]   = data_H[ncells-1];
          ndata_xnuc[npart*3+1] = data_HE[ncells-1];
          ndata_xnuc[npart*3+2] = 1. - data_H[ncells-1] - data_HE[ncells-1];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[npart*nspecies+k] = data_xnuc[k];
        }
        
        npart++;
        
      }
    }
    
    free( box );
  }
  
  memset( ndata_vel, 0, 3 * npart * sizeof(double) );
  
  printf( "Created %d particles.\n", npart );

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "pos", (PyObject*)createPyArray( ndata_pos, npart, 3 ) );
  free( ndata_pos );
  PyDict_SetStolenItem( dict, "mass", (PyObject*)createPyArray( ndata_mass, npart, 1 ) );
  free( ndata_mass );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( ndata_u, npart, 1 ) );
  free( ndata_u );
  PyDict_SetStolenItem( dict, "vel", (PyObject*)createPyArray( ndata_vel, npart, 3 ) );
  free( ndata_vel );
  PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( ndata_xnuc, npart, nspecies ) );
  free( ndata_xnuc );
  PyDict_SetStolenItem( dict, "count", (PyObject*)PyLong_FromLong( npart ) );
  
  return dict;
}

PyObject* _create_particles_healpix_shells(PyObject *self, PyObject *args, PyObject *keywds) {
  PyObject *data, *dict;
  t_helm_eos_table *helm_eos_table;
  PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
  double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
  double mass, mtot, vol, irho, iu, *ixnuc;
  int nspecies, nangular, ncells, nshells, nc, n, np, shell;
  double *ndata_pos, *ndata_mass, *ndata_u, *ndata_vel, *ndata_xnuc;
  int arraylen, i, j, k, p, index;
  double maxradius, radius, dr, vec[3];
  int makebox, boxres, nbox, *box, ix, iy, iz, idx;
  double boxsize, hboxsize, boxcellsize, cx, cy, cz, px, py, pz, minenergy;
  int writedensity, randomizeshells, noxnuc;
  double phi1, phi2, phi3, x, y, z, x2, y2, z2;

  static char *kwlist[] = { "data", "eos", "nshells", "nangular", "nspecies", "makebox", "boxsize", "boxres", "writedensity", "randomizeshells", "minenergy", NULL };
  
  nangular = 0;
  nspecies = 3;
  makebox = 0;
  boxsize = 0;
  boxres = 16;
  writedensity = 0;
  randomizeshells = 0;
  minenergy = 0;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "OO&i|iiidiiiid:create_particles_healpix_shells( data, eos, nshells, [nangular, nspecies, makebox, boxsize, boxres, writedensity, randomizeshells, minenergy] )", kwlist, &data, &pyConvertHelmEos, &helm_eos_table, &nshells, &nangular, &nspecies, &makebox, &boxsize, &boxres, &writedensity, &randomizeshells, &minenergy )) {
    return 0;
  }
  
  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  rho = (PyArrayObject*)PyDict_GetItemString( data, "rho" );
  data_rho = (double*)PyArray_DATA(rho);
  u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
  data_u = (double*)PyArray_DATA(u);

  data_H = data_HE = data_xnuc = NULL;
  if (nspecies == 3 && !PyDict_Contains( data, PyUnicode_FromString( "xnuc" ))) {
    noxnuc = 1;
    H = (PyArrayObject*)PyDict_GetItemString( data, "H" );
    data_H = (double*)PyArray_DATA(H);
    HE = (PyArrayObject*)PyDict_GetItemString( data, "HE" );
    data_HE = (double*)PyArray_DATA(HE);
  } else {
    noxnuc = 0;
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    data_xnuc = (double*)PyArray_DATA(xnuc);
  }

  seed();

  ixnuc = malloc( nspecies * sizeof(double) );

  mtot = 0.0;
  for (i=0; i<ncells; i++) mtot += data_dm[i];
  printf( "Total mass: %g solar masses\n", mtot / 1.989e33 );

  maxradius = data_r[ ncells-1 ];
  dr = maxradius / nshells;
  printf( "Using %d radial shells with a radial separation of %g.\n", nshells, dr );
  
  arraylen = 0x10000;
  ndata_pos = (double*)malloc( 3*arraylen*sizeof(double) );
  ndata_mass = (double*)malloc( arraylen*sizeof(double) );
  ndata_u = (double*)malloc( arraylen*sizeof(double) );
  ndata_vel = (double*)malloc( 3*arraylen*sizeof(double) );
  ndata_xnuc = (double*)malloc( nspecies*arraylen*sizeof(double) );
  p = 0;

  n = 0;
  if (nangular > 0) {
    /* nangunlar is the number of bins in phi-direction */
    nc = ( 0.5 * nangular ) * nangular;
    n = floor( sqrt( nc / 12. ) + 0.5 );
    printf( "Using %d particles per shell.\n", 12*n*n );
  }

  index = 0;
  mass = 0.;
  phi1 = phi2 = phi3 = 0;
  
  for (shell=0; shell<nshells; shell++) {
    radius = (shell+0.5) * dr;
    vol = 4./3. * M_PI * ( pow((shell+1.)*dr,3) - pow(shell*dr,3) );
    
    while (data_r[index] < radius)
      index++;
    irho = interpolate1D( index, radius, ncells, data_r, data_rho );
    iu = interpolate1D( index, radius, ncells, data_r, data_u );

    if (noxnuc) {
      ixnuc[0] = interpolate1D( index, radius, ncells, data_r, data_H );
      ixnuc[1] = interpolate1D( index, radius, ncells, data_r, data_HE );
      ixnuc[2] = 1. - ixnuc[0] - ixnuc[1];
    } else {
      for (i=0; i<nspecies; i++)
        ixnuc[i] = data_xnuc[i];
    }

    mass += irho * vol;

    if (nangular == 0) {
      nc = 4. * M_PI * radius * radius / (dr * dr);
      n = floor( sqrt( nc / 12. ) + 0.5 );
    }

    if (randomizeshells) {
      phi1 = randMT() * 2. * M_PI;
      phi2 = randMT() * 2. * M_PI;
      phi3 = randMT() * 2. * M_PI;
    }

    np = 12 * n * n;
    for (j=0; j<np; j++) {
      pix2vec_ring( n, j, vec );
      
      if (randomizeshells) {
        x = vec[0];
        y = cos( phi1 )*vec[1] - sin( phi1 )*vec[2];
        z = sin( phi1 )*vec[1] + cos( phi1 )*vec[2];
        
        x2 =  cos( phi2 )*x + sin( phi2 )*z;
        y2 = y;
        z2 = -sin( phi2 )*x + cos( phi2 )*z;
          
        ndata_pos[p*3+0] = radius*( cos( phi3 )*x2 - sin( phi3 )*y2 );
        ndata_pos[p*3+1] = radius*( sin( phi3 )*x2 + cos( phi3 )*y2 );
        ndata_pos[p*3+2] = radius*z2;
      } else {
        for (k=0; k<3; k++)
    ndata_pos[p*3+k] = vec[k] * radius;
      }
      
      if (writedensity)
        ndata_mass[p] = irho;
      else
        ndata_mass[p] = irho * vol / np;
      ndata_u[p] = iu;

      for (k=0; k<nspecies; k++)
        ndata_xnuc[p*nspecies+k] = ixnuc[k];
      
      p++;

      if (p == arraylen) {
        ndata_pos = resize( ndata_pos, arraylen*3, 2*arraylen*3 );
        ndata_vel = resize( ndata_vel, arraylen*3, 2*arraylen*3 );
        ndata_mass = resize( ndata_mass, arraylen, 2*arraylen );
        ndata_u = resize( ndata_u, arraylen, 2*arraylen );
        ndata_xnuc = resize( ndata_xnuc, arraylen*nspecies, 2*arraylen*nspecies );
        arraylen *= 2;
      }
    }   
  }
  
  if (makebox) {
    nbox = boxres * boxres * boxres;
    box = (int*)malloc( nbox * sizeof(int) );
    memset( box, 0, nbox * sizeof(int) );
    
    hboxsize = boxsize / 2.;
    boxcellsize = boxsize / (double)boxres;
    for (i=0; i<p; i++) {
      ndata_pos[i*3+0] += hboxsize;
      ndata_pos[i*3+1] += hboxsize;
      ndata_pos[i*3+2] += hboxsize;
      
      px = ndata_pos[i*3+0];
      py = ndata_pos[i*3+1];
      pz = ndata_pos[i*3+2];
      
      ix = min( max( floor( px / boxcellsize ), 0 ), boxres-1 );
      iy = min( max( floor( py / boxcellsize ), 0 ), boxres-1 );
      iz = min( max( floor( pz / boxcellsize ), 0 ), boxres-1 );
      
      cx = (ix+0.5) * boxcellsize;
      cy = (iy+0.5) * boxcellsize;
      cz = (iz+0.5) * boxcellsize;

      idx = 0;
      if (px > cx) idx += 1;
      if (py > cy) idx += 2;
      if (pz > cz) idx += 4;

      box[(boxres*iz+iy)*boxres+ix] |= 1 << idx;
    }
    
    for (i=0; i<nbox; i++) {
      if (box[i] != 0xFF) {
        ix = i % boxres;
        iy = (i / boxres) % boxres;
        iz = i / (boxres*boxres);
        
        ndata_pos[p*3+0] = (ix+0.5) * boxcellsize;
        ndata_pos[p*3+1] = (iy+0.5) * boxcellsize;
        ndata_pos[p*3+2] = (iz+0.5) * boxcellsize;
        
        ndata_mass[p] = 0;
        ndata_u[p] = minenergy;
        if (noxnuc) {
          ndata_xnuc[p*3]   = data_H[ncells-1];
          ndata_xnuc[p*3+1] = data_HE[ncells-1];
          ndata_xnuc[p*3+2] = 1. - data_H[ncells-1] - data_HE[ncells-1];
        } else {
          for (k=0; k<nspecies; k++)
            ndata_xnuc[p*nspecies+k] = data_xnuc[k];
        }
        
        p++;
        
        if (p == arraylen) {
          ndata_pos = resize( ndata_pos, arraylen*3, 2*arraylen*3 );
          ndata_vel = resize( ndata_vel, arraylen*3, 2*arraylen*3 );
          ndata_mass = resize( ndata_mass, arraylen, 2*arraylen );
          ndata_u = resize( ndata_u, arraylen, 2*arraylen );
          ndata_xnuc = resize( ndata_xnuc, arraylen*nspecies, 2*arraylen*nspecies );
          arraylen *= 2;
        }
      }
    }
    
    free( box );
  }

  free( ixnuc );
  
  memset( ndata_vel, 0, 3 * p * sizeof(double) );
  
  printf( "Created %d particles.\n", p );

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "pos", (PyObject*)createPyArray( ndata_pos, p, 3 ) );
  free( ndata_pos );
  PyDict_SetStolenItem( dict, "mass", (PyObject*)createPyArray( ndata_mass, p, 1 ) );
  free( ndata_mass );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( ndata_u, p, 1 ) );
  free( ndata_u );
  PyDict_SetStolenItem( dict, "vel", (PyObject*)createPyArray( ndata_vel, p, 3 ) );
  free( ndata_vel );
  PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( ndata_xnuc, p, nspecies ) );
  free( ndata_xnuc );
  PyDict_SetStolenItem( dict, "count", (PyObject*)PyLong_FromLong( p ) );
  
  return dict;
}

PyObject* _create_particles_random_shells(PyObject *self, PyObject *args) {
  FILE *fd;
  PyObject *data;
  char *filename, fn[200];
  double pmass, mass;
  int convert2HHeC, ncells, nspecies, nspecies_orig, ncell, nparttot, npart, npartmax, npartcell, npartcellmax, p, ptot, pmax, filenr, usecells, cellspercell, index, speciesdim;
  PyArrayObject *dm, *r, *v, *u, *xnuc;
  double *data_dm, *data_r, *data_v, *data_u, *data_xnuc;
  int i, split, num_files;
  float *ndata_pos, *ndata_vel, *ndata_mass, *ndata_u, *ndata_xnuc;
  int *ndata_id;
  double last_radius, next_radius, radius, phi, theta, vel, ene, xHe;

  npartmax = MAXPARTPERICFILE; /* maximum number of particles in on file. if number of particles is larger, file is split */

  usecells = 0;
  convert2HHeC = 0;
  if (!PyArg_ParseTuple( args, "Ods|ii:create_particles_random_shells( data, pmass, filename, [usecells, convert2HHeC] )", &data, &pmass, &filename, &usecells, &convert2HHeC )) {
    Py_RETURN_FALSE;
  }

  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  ncells = PyArray_DIMS(dm)[0];
  data_dm = (double*)PyArray_DATA(dm);
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  v = (PyArrayObject*)PyDict_GetItemString( data, "v" );
  data_v = (double*)PyArray_DATA(v);

  data_u = data_xnuc = NULL;
  if (PyDict_Contains( data, PyUnicode_FromString( "u" ))) {
    u = (PyArrayObject*)PyDict_GetItemString( data, "u" );
    data_u = (double*)PyArray_DATA(u);
  }

  speciesdim = 1;
  if (PyDict_Contains( data, PyUnicode_FromString( "xnuc" ))) {
    xnuc = (PyArrayObject*)PyDict_GetItemString( data, "xnuc" );
    nspecies_orig = PyArray_DIMS(xnuc)[0];
    if (PyArray_NDIM(xnuc) > 1) {
      speciesdim = PyArray_DIMS(xnuc)[1];
      if (speciesdim != ncells) {
        printf( "xnuc has to be of dimension [nspecies] or [nspecies,ncells]\n" );
        Py_RETURN_FALSE;
      }
    }
    data_xnuc = (double*)PyArray_DATA(xnuc);
  } else {
    nspecies_orig = 0;
    speciesdim = 1;
  }

  if (convert2HHeC) {
    nspecies = 3;
  } else {
    nspecies = nspecies_orig;
  }
  printf( "Particles have %d species.\n", nspecies );

  if (!usecells || usecells > ncells) {
    usecells = ncells;
  }

  cellspercell = max( 1, floor( ncells / usecells + 0.5 ) );
  if (usecells && cellspercell > 1) {
    usecells = floor( ncells / cellspercell );
    printf( "Cells per cell: %d, using %d artificial cells.\n", cellspercell, usecells );
  }

  mass = 0;
  for (i=0; i<ncells; i++) {
    mass += data_dm[i];
  }
  nparttot = floor( mass / pmass + 0.5 );

  split = nparttot > npartmax ? 1 : 0;
  num_files = ceil( nparttot / npartmax );

  npart = min( nparttot, npartmax );
  ndata_pos = (float*)malloc( 3.0 * npart * sizeof(float) );
  ndata_vel = (float*)malloc( 3.0 * npart * sizeof(float) );
  ndata_mass = (float*)malloc( npart * sizeof(float) );
  ndata_u = (float*)malloc( npart * sizeof(float) );
  ndata_id = (int*)malloc( npart * sizeof(int) );
  ndata_xnuc = (float*)malloc( nspecies * npart * sizeof(float) );
  
  seed(); /* initialize random number generator */

  ptot = 0;
  pmax = 0;
  filenr = 0;
  ncell = -1;
  npartcellmax = 0;
  next_radius = 0;
  last_radius = 0;
  mass = 0;

  do {
    npart = min( max( 0, nparttot-npartmax*filenr ), npartmax );
    if (!npart) {
      break;
    }

    if (split) {
      sprintf( fn, "%s.%03d", filename, filenr );
    } else {
      sprintf( fn, "%s", filename );
    }
    p = 0;
    pmax += npart;

    while (ncell < usecells && ptot < pmax) {
      while (ptot < npartcellmax && ptot < pmax) {
        radius = pow( randMT() * ( pow(next_radius,3.) - pow(last_radius,3.) ) + pow(last_radius,3.), 1./3. );
        phi = randMT() * 2. * M_PI;
        theta = acos( 2. * randMT() - 1. );

        index = ncell * cellspercell;
        while (data_r[index] < radius)
          index++;

        if (index >= ncells || index < 1) printf( "%d: %g %g %g\n", index, last_radius, radius, next_radius );

        ndata_pos[p*3]   = radius * sin( theta ) * cos( phi );
        ndata_pos[p*3+1] = radius * sin( theta ) * sin( phi );
        ndata_pos[p*3+2] = radius * cos( theta );
        
        ndata_mass[p] = pmass;
        ndata_id[p] = ptot+1;
        
        vel = interpolate1D( index, radius, ncells, data_r, data_v );
        ndata_vel[p*3]   = vel * ndata_pos[p*3]   / radius;
        ndata_vel[p*3+1] = vel * ndata_pos[p*3+1] / radius;
        ndata_vel[p*3+2] = vel * ndata_pos[p*3+2] / radius;

        if (data_u) {
          ene = interpolate1D( index, radius, ncells, data_r, data_u );
          ndata_u[p] = max( ene - 0.5 * vel*vel, 0.001 * 0.5 * vel*vel );
        } else {
          ndata_u[p] = 0.001 * 0.5 * vel*vel;
        }

        if (nspecies) {
          if (convert2HHeC) {
            if (nspecies_orig) {
              xHe = interpolate1D_valuedim( index, radius, ncells, data_r, data_xnuc, nspecies_orig, 0 );
            } else {
              xHe = 0;
            }
            ndata_xnuc[p*3]   = 0;
            ndata_xnuc[p*3+1] = xHe;
            ndata_xnuc[p*3+2] = 1.0 - xHe;
          } else {
            if (speciesdim > 1) {
              for (i=0; i<nspecies; i++) {
                ndata_xnuc[p*nspecies+i] = interpolate1D_valuedim( index, radius, ncells, data_r, data_xnuc, nspecies, i );
              }
            } else {
              for (i=0; i<nspecies; i++) {
                ndata_xnuc[p*nspecies+i] = data_xnuc[i];
              }
            }
          }
        }

        p++;
        ptot++;
        mass -= pmass;
      }
      
      if (ptot == npartcellmax) {
        ncell++;

        for (i=ncell*cellspercell; i<(ncell+1)*cellspercell; i++)
          mass += data_dm[i];

        npartcell = floor( mass / pmass + 0.5 );
        npartcellmax += npartcell;

        last_radius = next_radius;
        next_radius = data_r[ (ncell+1)*cellspercell ];
      }
    }

    fd = fopen( fn, "w" );

    gadget_writeHeader( fd, nparttot, npart, num_files );
    gadget_writeBlockFloat( fd, "POS ", npart, 3, ndata_pos );
    gadget_writeBlockFloat( fd, "VEL ", npart, 3, ndata_vel );
    gadget_writeBlockInt(   fd, "ID  ", npart, 1, ndata_id );
    gadget_writeBlockFloat( fd, "MASS", npart, 1, ndata_mass );
    gadget_writeBlockFloat( fd, "U   ", npart, 1, ndata_u );
    gadget_writeBlockFloat( fd, "XNUC", npart, nspecies, ndata_xnuc );

    fclose( fd );

    filenr++;
  } while (1);

  free( ndata_pos );
  free( ndata_vel );
  free( ndata_mass );
  free( ndata_u );
  free( ndata_id );
  free( ndata_xnuc );

  Py_RETURN_TRUE;
}

PyObject* _create_particles_random_tori(PyObject *self, PyObject *args) {
  FILE *fd;
  PyObject *data;
  char *filename, fn[200];
  double pmass, rmin;
  int mirror, nspecies, nr, nz, p, ptot, filenr, npartmax, reverse;
  PyArrayObject *dm, *r, *z, *vr, *vz;
  double *data_dm, *data_r, *data_z, *data_vr, *data_vz;
  int ir, iz, in, split, num_files, nparttot, npart;
  double r1, r2, z1, z2;
  float *ndata_pos, *ndata_vel, *ndata_mass, *ndata_u, *ndata_xnuc;
  int *ndata_id;
  double rr, zz, phi, velr, velz;

  npartmax = MAXPARTPERICFILE; /* maximum number of particles in on file. if number of particles is larger, file is split */

  mirror = 0;
  rmin = 0;
  if (!PyArg_ParseTuple( args, "Ods|id:create_particles_random_tori( data, pmass, filename, [mirror, rmin] )", &data, &pmass, &filename, &mirror, &rmin )) {
    Py_RETURN_FALSE;
  }

  dm = (PyArrayObject*)PyDict_GetItemString( data, "dm" );
  nr = PyArray_DIMS(dm)[0];
  nz = PyArray_DIMS(dm)[1];
  if (PyArray_STRIDES(dm)[0] == 8 && PyArray_STRIDES(dm)[1] == nr*8) {
    reverse = 0;
  } else if (PyArray_STRIDES(dm)[0] == nz*8 && PyArray_STRIDES(dm)[1] == 8) {
    reverse = 1;
  } else {
    reverse = 0;
    printf( "Corrupt array data[\"dm\"], stopping.\n" );
  }

  data_dm = (double*)PyArray_DATA(dm);

  /* r and z are assumed to have a length of nr+1 and nz+1 */
  r = (PyArrayObject*)PyDict_GetItemString( data, "r" );
  data_r = (double*)PyArray_DATA(r);
  z = (PyArrayObject*)PyDict_GetItemString( data, "z" );
  data_z = (double*)PyArray_DATA(z);
  
  vr = (PyArrayObject*)PyDict_GetItemString( data, "vr" );
  data_vr = (double*)PyArray_DATA(vr);
  vz = (PyArrayObject*)PyDict_GetItemString( data, "vz" );
  data_vz = (double*)PyArray_DATA(vz);

  nspecies = 3;

  nparttot = 0;
  for (ir=0; ir<nr; ir++) {
                r1 = data_r[ir];
    for (iz=0; iz<nz; iz++) {
            z1 = data_z[iz];
      if (sqrt(r1*r1+z1*z1) >= rmin)
         nparttot += floor( data_dm[iz*nr+ir] / pmass + 0.5 );
    }
  }

  if (mirror) {
    nparttot *= 2;
  }

  split = nparttot > npartmax ? 1 : 0;
  num_files = ceil( (double)nparttot / (double)npartmax );

  printf( "Writing ICs to %d files.\n", num_files );

  npart = min( nparttot, npartmax );
  ndata_pos = (float*)malloc( 3.0 * npart * sizeof(float) );
  ndata_vel = (float*)malloc( 3.0 * npart * sizeof(float) );
  ndata_mass = (float*)malloc( npart * sizeof(float) );
  ndata_u = (float*)malloc( npart * sizeof(float) );
  ndata_id = (int*)malloc( npart * sizeof(int) );
  ndata_xnuc = (float*)malloc( nspecies * npart * sizeof(float) );
  
  seed(); /* initialize random number generator */

  p = 0;
  ptot = 0;
  filenr = 0;

  for (ir=0; ir<nr; ir++) {
    r1 = data_r[ir];
    r2 = data_r[ir+1];
    for (iz=0; iz<nz; iz++) {
      z1 = data_z[iz];
      z2 = data_z[iz+1];

      if (sqrt(r1*r1+z1*z1)< rmin)
        continue;

      if (!reverse) {
        npart = floor( data_dm[iz*nr+ir] / pmass + 0.5 );
      } else {
        npart = floor( data_dm[ir*nz+iz] / pmass + 0.5 );
      }
      for (in=0; in<npart; in++) {
        if (p == npartmax || (mirror && p+1 == npartmax)) {
          /* we always need more than one file, if we reach this point */
          sprintf( fn, "%s.%d", filename, filenr );
          
          fd = fopen( fn, "w" );

          gadget_writeHeader( fd, nparttot, p, num_files );
          gadget_writeBlockFloat( fd, "POS ", p, 3, ndata_pos );
          gadget_writeBlockFloat( fd, "VEL ", p, 3, ndata_vel );
          gadget_writeBlockInt(   fd, "ID  ", p, 1, ndata_id );
          gadget_writeBlockFloat( fd, "MASS", p, 1, ndata_mass );
          gadget_writeBlockFloat( fd, "U   ", p, 1, ndata_u );
          gadget_writeBlockFloat( fd, "XNUC", p, nspecies, ndata_xnuc );

          fclose( fd );

          filenr++;
          p = 0;
        }

        rr = pow( randMT() * ( r2*r2 - r1*r1 ) + r1*r1, 0.5 );
        phi = randMT() * 2. * M_PI;
        zz = randMT() * (z2 - z1) + z1;

        ndata_pos[p*3]   = rr * cos( phi );
        ndata_pos[p*3+1] = rr * sin( phi );
        ndata_pos[p*3+2] = zz;

        ndata_mass[p] = pmass;
        ndata_id[p] = ptot+1; /* IDs start with 1 */

        velr = interpolate2D( rr, zz, data_r, data_z, data_vr, nr, nz, reverse );
        velz = interpolate2D( rr, zz, data_r, data_z, data_vz, nr, nz, reverse );

        ndata_vel[p*3]   = velr * ndata_pos[p*3]   / rr;
        ndata_vel[p*3+1] = velr * ndata_pos[p*3+1] / rr;
        ndata_vel[p*3+2] = velz;

        ndata_u[p] = 0.001 * 0.5 * (velr*velr + velz*velz);

        ndata_xnuc[p*3]   = 0;
        ndata_xnuc[p*3+1] = 0;
        ndata_xnuc[p*3+2] = 1;

        p++;
        ptot++;

        if (!mirror)
          continue;

        /* otherwise do mirror particle */
        rr = pow( randMT() * ( r2*r2 - r1*r1 ) + r1*r1, 0.5 );
        phi = randMT() * 2. * M_PI;
        zz = randMT() * (z2 - z1) + z1;

        ndata_pos[p*3]   = rr * cos( phi );
        ndata_pos[p*3+1] = rr * sin( phi );
        ndata_pos[p*3+2] = -zz;

        ndata_mass[p] = pmass;
        ndata_id[p] = ptot+1; /* IDs start with 1 */

        velr = interpolate2D( rr, zz, data_r, data_z, data_vr, nr, nz, reverse );
        velz = interpolate2D( rr, zz, data_r, data_z, data_vz, nr, nz, reverse );

        ndata_vel[p*3]   = velr * ndata_pos[p*3]   / rr;
        ndata_vel[p*3+1] = velr * ndata_pos[p*3+1] / rr;
        ndata_vel[p*3+2] = -velz;

        ndata_u[p] = 0.001 * 0.5 * (velr*velr + velz*velz);

        ndata_xnuc[p*3]   = 0;
        ndata_xnuc[p*3+1] = 0;
        ndata_xnuc[p*3+2] = 1;

        p++;
        ptot++;
      }
    }
  }

  if (p > 0) {
    if (split) {
      sprintf( fn, "%s.%d", filename, filenr );
    } else {
      sprintf( fn, "%s", filename );
    }
  
    fd = fopen( fn, "w" );

    gadget_writeHeader( fd, nparttot, p, num_files );
    gadget_writeBlockFloat( fd, "POS ", p, 3, ndata_pos );
    gadget_writeBlockFloat( fd, "VEL ", p, 3, ndata_vel );
    gadget_writeBlockInt(   fd, "ID  ", p, 1, ndata_id );
    gadget_writeBlockFloat( fd, "MASS", p, 1, ndata_mass );
    gadget_writeBlockFloat( fd, "U   ", p, 1, ndata_u );
    gadget_writeBlockFloat( fd, "XNUC", p, nspecies, ndata_xnuc );

    fclose( fd );
  }

  free( ndata_pos );
  free( ndata_vel );
  free( ndata_mass );
  free( ndata_u );
  free( ndata_id );
  free( ndata_xnuc );

  Py_RETURN_TRUE;
}


PyObject* _create_particles_fill_grid(PyObject *self, PyObject *args, PyObject *keywds) {
  PyArrayObject *pyPos;
  double *pos, *grid_new;
  double boxsize, cellsize;
  double boxx, boxy, boxz;
  int *grid;
  int npart, part, pp, ncells, boxres;
  int ix, iy, iz, idx, nx, ny, nz;
  double px, py, pz;
  double longx, longy, longz;

  static char *kwlist[] = { "pos", "boxsize", "boxres", "longx", "longy", "longz", NULL };
  
  longx = 1.0;
  longy = 1.0;
  longz = 1.0;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "O!di|ddd:create_particles_fill_grid( pos, boxsize, boxres, [longx, longy, longz] )", kwlist, &PyArray_Type, &pyPos, &boxsize, &boxres, &longx, &longy, &longz )) {
    return 0;
  }

  if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [:,3] and type double" );
    return 0;
        }

  pos = (double*)PyArray_DATA(pyPos);
  npart = PyArray_DIMS(pyPos)[0];

  printf( "%d particles found.\n", npart );

  cellsize = boxsize / boxres;
  boxx = boxsize * longx;
  boxy = boxsize * longy;
  boxz = boxsize * longz;

  nx = floor(boxx / cellsize);
  ny = floor(boxy / cellsize);
  nz = floor(boxz / cellsize);

  printf( "Building grid with %d x %d x %d cells.\n", nx, ny, nz );
  ncells = nx * ny * nz;

  grid = malloc( ncells * sizeof(int) );
  memset( grid, 0, ncells*sizeof(int) );
  
  for (part=0; part<npart; part++) {
    px = *(double*)((char*)pos + part*PyArray_STRIDES(pyPos)[0] + 0*PyArray_STRIDES(pyPos)[1]);
    py = *(double*)((char*)pos + part*PyArray_STRIDES(pyPos)[0] + 1*PyArray_STRIDES(pyPos)[1]);
    pz = *(double*)((char*)pos + part*PyArray_STRIDES(pyPos)[0] + 2*PyArray_STRIDES(pyPos)[1]);

    ix = floor( px / cellsize );
    iy = floor( py / cellsize );
    iz = floor( pz / cellsize );

    if (ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz) {
      grid[ ((iz * ny + iy) * nx) + ix ] += 1;
    }
  }

  pp = 0;
  grid_new = malloc( ncells * sizeof(double) * 3 );

  idx = 0;
  
  for (ix=0; ix<nx; ix++)
    for (iy=0; iy<ny; iy++)
      for (iz=0; iz<nz; iz++) {
        idx = ((iz * ny + iy) * nx) + ix;
        if (grid[idx] == 0) {
          grid_new[ pp*3   ] = (ix+0.5) * cellsize;
          grid_new[ pp*3+1 ] = (iy+0.5) * cellsize;
          grid_new[ pp*3+2 ] = (iz+0.5) * cellsize;

          pp++;
        }
      }

  free( grid );
  grid_new = realloc( grid_new, pp * sizeof(double) * 3 );
  
  printf( "Created %d particles to fill grid.\n", pp );
  
        PyObject* res = (PyObject*)createPyArray( grid_new, pp, 3 );
  free( grid_new );

  return res;
}

PyObject* _create_tracers(PyObject *self, PyObject *args, PyObject *keywds) {
  PyArrayObject *pyRho;
  double *rho, rhomax, irhomax;
  long ntracer, ndone;
  int nx, ny, nz;
  long ix, iy, iz, i, ncells;
  double px, py, pz, *pos;
  int seed, nextoutput, slowly = 0;
  double runtime;
  clock_t start;

  static char *kwlist[] = { "rho", "ntracer", "seed", NULL };

  start = clock();

  seed = 4357;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "O!l|i:createTracers( rho, ntracer, [seed] )", kwlist, &PyArray_Type, &pyRho, &ntracer, &seed )) {
    return 0;
  }

  if (PyArray_NDIM(pyRho) != 3 || PyArray_TYPE(pyRho) != NPY_DOUBLE) {
    PyErr_SetString( PyExc_ValueError, "rho has to be of dimensions [:,:,:] and type double" );
    return 0;
  }

  seedMT( seed );

  rho = (double*)PyArray_DATA(pyRho);
  nx = PyArray_DIMS(pyRho)[0];
  ny = PyArray_DIMS(pyRho)[1];
  nz = PyArray_DIMS(pyRho)[2];
  ncells = nx * ny * nz;

  rhomax = 0;
  for (i = 0; i < ncells; i++)
    if (rho[i] > rhomax)
      rhomax = rho[i];

  irhomax = 1. / rhomax;
  for (i = 0; i < ncells; i++)
    rho[i] *= irhomax;

  pos = (double*)malloc( 3 * ntracer * sizeof(double) );

  printf( "Creating %ld tracers...", ntracer );

  ndone = 0;
  nextoutput = ntracer / 100;
  while (ndone < ntracer)
    {
      px = randMT();
      py = randMT();
      pz = randMT();

      ix = min( floor( px * nx ), nx-1 );
      iy = min( floor( py * ny ), ny-1 );
      iz = min( floor( pz * nz ), nz-1 );

      if ( randMT() < rho[ (ix * ny + iy) * nz + iz ] )
        {
          pos[ndone*3  ] = px;
          pos[ndone*3+1] = py;
          pos[ndone*3+2] = pz;

          ndone++;

          if (ndone >= nextoutput) {
            runtime = ((double)clock()-(double)start)/CLOCKS_PER_SEC;
                                  
            if (nextoutput == ntracer / 100) {
              if ( runtime > 60. ) {
                slowly = 1;
              } else {
                nextoutput = 0;
              }
            }

            printf( "%ld / %ld tracer done (%d%%): %ds elapsed, ~%ds remaining\n", ndone, ntracer, (int)floor(100.0*(double)ndone/(double)ntracer), (int)(runtime), (int)(runtime/ndone*(ntracer-ndone)) );

            if (slowly)
              nextoutput += ntracer / 100;
            else
              nextoutput += ntracer /  10;
          }
        }
    }

  printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );

  PyObject* res = (PyObject*)createPyArray( pos, ntracer, 3 );
  free( pos );

  return res;
}

PyObject* _healpix2vec_ring(PyObject *self, PyObject *args ) {
  int n, j;
  double vec[3];

  if (!PyArg_ParseTuple( args, "ii:healpix2vec_ring( n, j )", &n, &j )) {
    return 0;
  }

  pix2vec_ring( n, j, vec );
  return (PyObject*)createPyArray( vec, 3, 1 );
}


static PyMethodDef createICsMethods[] = {
  { "create_particles_cube", _create_particles_cube, METH_VARARGS, "" },
  { "create_particles_healpix", (PyCFunction)_create_particles_healpix, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_particles_healpix_vol", (PyCFunction)_create_particles_healpix_vol, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_particles_healpix_grad", (PyCFunction)_create_particles_healpix_grad, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_particles_healpix_shells", (PyCFunction)_create_particles_healpix_shells, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_particles_random_shells", _create_particles_random_shells, METH_VARARGS, "" },
  { "create_particles_random_tori", _create_particles_random_tori, METH_VARARGS, "" },
  { "create_particles_fill_grid", (PyCFunction)_create_particles_fill_grid, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_tracers", (PyCFunction)_create_tracers, METH_VARARGS|METH_KEYWORDS, "" },
  { "healpix2vec_ring", _healpix2vec_ring, METH_VARARGS, "" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "createICs", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  createICsMethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_createICs(void)
{
  import_array();
  return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initcreateICs(void)
{
  Py_InitModule( "createICs", createICsMethods );
  import_array();
}
#endif

