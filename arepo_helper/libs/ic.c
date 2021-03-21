#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helm_eos.h"
#include "pyhelm_eos.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#define G 6.6738e-8 /* gravitational constant */

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

inline void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
  PyDict_SetItemString(dict, key, object);
  Py_DECREF(object);
}

PyArrayObject* createPyArray( double *data, int length ) {
  PyArrayObject* pyData;
  
  npy_intp dims[1];
  dims[0] = length;
  
  pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_DOUBLE );
  memcpy( PyArray_DATA(pyData), data, length*sizeof(double) );
  
  return pyData;
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

struct paramsWD {
  double temp;
  double *xnuc;
  t_helm_eos_table *eos;
  double rho;
};

int createWDIntegrator( double r, const double y[], double f[], void *params)
{
  double temp = ((struct paramsWD*)params)->temp;
  double *xnuc = ((struct paramsWD*)params)->xnuc;
  t_helm_eos_table *eos = ((struct paramsWD*)params)->eos;
  double rho = ((struct paramsWD*)params)->rho;
  
  struct eos_result res;
  eos_calc_ptgiven( eos, y[0], xnuc, temp, &rho, &res );  
  r = max( r, 1e2 );
  double r2 = r*r;
  
  f[0] = - G * y[1] * rho / r2;
  f[1] = 4. * M_PI  * r2 * rho;
  
  ((struct paramsWD*)params)->rho = rho;
  
  return GSL_SUCCESS;
}

int createWDIntegratorShell( double r, const double h[], double f[], void *params)
{
  //double temp = ((struct paramsWD*)params)->temp;
  double temp = h[2];
  double *xnuc = ((struct paramsWD*)params)->xnuc;
  t_helm_eos_table *eos = ((struct paramsWD*)params)->eos;
  double rho = ((struct paramsWD*)params)->rho;

  struct eos_result res;
  eos_calc_ptgiven( eos, h[0], xnuc, temp, &rho, &res );
  r = max( r, 1e2 );
  double r2 = r*r;

  f[0] = - G * h[1] * rho / r2;
  f[1] = 4. * M_PI  * r2 * rho;
  f[2] = f[0] * temp / h[0] * (res.gamma_2 - 1) / res.gamma_2;
//  printf("%e %e %e\n", f[2], h[2], ((struct paramsWD*)params)->temp);

  ((struct paramsWD*)params)->rho = rho;

  return GSL_SUCCESS;
}

PyObject* _createWhiteDwarf(PyObject *self, PyObject *args, PyObject *keywds) {
  double rho0, temp; /* central density and temperature */
  double xHe4, xC12, xO16, xNe20, xMg24; /* mass fractions */
  double tol;
  t_helm_eos_table *eos;
  double *r, *p, *e, *rho, *dm, *csnd;
  int i, arraylen, count;
  PyObject* dict;
  struct eos_result res;
  
  static char *kwlist[] = { "eos", "rho0", "temp", "xHe4", "xC12", "xO16", "xNe20", "xMg24", "tol", NULL };
  
  tol = 1e-6;
  temp = 5e5;
  xHe4 = xC12 = xO16 = xNe20 = xMg24 = 0.;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "O&d|ddddddd:create_wd( eos, rho0, [temp, xHe4, xC12, xO16, xNe20, xMg24, tol] )", kwlist, &pyConvertHelmEos, &eos, &rho0, &temp, 
    &xHe4, &xC12, &xO16, &xNe20, &xMg24, &tol )) {
    return 0;
  }
  
  arraylen = 0x10000;
  r = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  rho = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );
  csnd = (double*)malloc( arraylen * sizeof(double) );

  double xtot = xHe4 + xC12 + xO16 + xNe20 + xMg24;
  printf("Abundances: He4=%g, C12=%g O16=%g Ne20=%g, Mg24=%g, sum=%g.\n", xHe4, xC12, xO16, xNe20, xMg24, xtot );
  if(fabs(xtot-1.0) > 1e-14)
    {
      PyErr_SetString( PyExc_ValueError, "Inconsistent Abundances.\n" );
      return 0;
    }
  
  double xnuc[eos->nspecies];
  for(i = 0; i < eos->nspecies; i++)
    {
      xnuc[i] = 0.;
      if(eos->nz[i] == 2 && eos->na[i] == 4) xnuc[i] = xHe4;
      if(eos->nz[i] == 6 && eos->na[i] == 12) xnuc[i] = xC12;
      if(eos->nz[i] == 8 && eos->na[i] == 16) xnuc[i] = xO16;
      if(eos->nz[i] == 10 && eos->na[i] == 20) xnuc[i] = xNe20;
      if(eos->nz[i] == 12 && eos->na[i] == 24) xnuc[i] = xMg24;
    }

  eos_calc_tgiven( eos, rho0, xnuc, temp, &res );
  rho[0] = rho0;
  r[0] = 1e2;
  p[0] = res.p.v;
  e[0] = res.e.v;
  dm[0] = rho0 * 4. / 3. * M_PI * r[0] * r[0] * r[0];
  csnd[0] = res.sound;
  count = 1;
  
  double y[2];
  y[0] = p[0];
  y[1] = 0.;
  
  struct paramsWD params;
  params.temp = temp;
  params.xnuc = xnuc;
  params.eos = eos;
  params.rho = rho0;
  
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, tol);
  gsl_odeiv_evolve * ev = gsl_odeiv_evolve_alloc (2);
  gsl_odeiv_system sys = {createWDIntegrator, 0, 2, &params};
  
  double mass = dm[0];
  double rad = 1e2;
  double drad = 1e3;
  while (rho[count-1] > 1e-5 && rad < 1e10)
    {
      int status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e10, &drad, y);
      
      if (status != GSL_SUCCESS)
        {
          PyErr_SetString( PyExc_ValueError, "ODE Solver failed.\n" );
          return 0;
        }
      
      /* increase arraylen if necessary */
      if (count == arraylen) {
        rho = resize( rho, arraylen, 2*arraylen );
        r = resize( r, arraylen, 2*arraylen );
        p = resize( p, arraylen, 2*arraylen );
        e = resize( e, arraylen, 2*arraylen );
        dm = resize( dm, arraylen, 2*arraylen );
        csnd = resize( csnd, arraylen, 2*arraylen );
        arraylen *= 2;
      }
      
      if(y[0] <= 0.)
        break;
      
      rho[count] = rho[count-1];
      if(rho[count] < 1e-5)
        break;
      
      eos_calc_ptgiven( eos, y[0], xnuc, temp, &rho[count], &res );
      r[count] = rad;
      p[count] = y[0];
      e[count] = res.e.v;
      dm[count] = y[1] - mass;
      csnd[count] = res.sound;
      mass = y[1];
      count++;
    }
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, count ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, count ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, count ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, count ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, count ) );
  PyDict_SetStolenItem( dict, "csnd", (PyObject*)createPyArray( csnd, count ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( count ) );
  
  free( r );
  free( p );
  free( e );
  free( rho );
  free( dm );
  free( csnd );

  gsl_odeiv_evolve_free(ev);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return dict;
}

PyObject* _createWhiteDwarfHeShell(PyObject *self, PyObject *args, PyObject *keywds) {
  double rho0, tempC, rhoS, tempS; /* central density and temperature, density at base of shell, temperature in shell */
  double xNe20, xMg24; /* mass fractions */
  double xHeC, xHeS, xCC, xCS, xOC, xOS; /* initial mass fractions for core and shell */
  t_helm_eos_table *eos;
  double *r, *p, *e, *rho, *dm;
  double *temp, *He4, *C12, *O16;
  int i, arraylen, count, var, trans;
  PyObject* dict;
  struct eos_result res;

  static char *kwlist[] = { "eos", "rho0", "tempC", "tempS", "rhoS", "xHeS", "xCS", "xOS", "xHeC", "xCC", "xOC", "xNe20", "xMg24", NULL};

  tempC = 5e5;
    tempS = 4e8; /* critical temperature for detonation: 4e9, get value for tempS from Fink 2010 */
  xHeC = xHeS = xCC = xCS = xOC = xOS = xNe20 = xMg24 = 0.;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "O&dddd|dddddddd:create_wdHeShell( eos, rho0, tempC, tempS, rhoS, xHeS, xCS, xOS, [xHeC, xCC, xOC, xNe20, xMg24] )", kwlist, &pyConvertHelmEos, &eos, &rho0, &tempC, &tempS, &rhoS, &xHeS, &xCS, &xOS, &xHeC, &xCC, &xOC, &xNe20, &xMg24)) {
    return 0;
  }

  arraylen = 0x10000;
  r = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  rho = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );
  temp = (double*)malloc( arraylen * sizeof(double) );
  He4 = (double*)malloc( arraylen * sizeof(double) );
  C12 = (double*)malloc( arraylen * sizeof(double) );
  O16 = (double*)malloc( arraylen * sizeof(double) );


  double xtotC = xHeC + xCC + xOC + xNe20 + xMg24;
  printf("Abundances in core: He4=%g, C12=%g O16=%g Ne20=%g, Mg24=%g, sum=%g.\n", xHeC, xCC, xOC, xNe20, xMg24, xtotC );
  if(fabs(xtotC-1.0) > 1e-14)
    {
      PyErr_SetString( PyExc_ValueError, "Inconsistent Abundances.\n" );
      return 0;
    }

  double xtotS = xHeS + xCS + xOS + xNe20 + xMg24;
  printf("Abundances in shell: He4=%g, C12=%g O16=%g Ne20=%g, Mg24=%g, sum=%g.\n", xHeS, xCS, xOS, xNe20, xMg24, xtotS );
  if(fabs(xtotS-1.0) > 1e-14)
    {
    PyErr_SetString( PyExc_ValueError, "Inconsistent Abundances.\n" );
    return 0;
    }

  double xHe4, xC12, xO16;
  xHe4 = xC12 = xO16 = 0.;
  xHe4 = xHeC;
  xC12 = xCC;
  xO16 = xOC;

  double xnuc[eos->nspecies];
  for(i = 0; i < eos->nspecies; i++)
    {
    xnuc[i] = 0.;
    if(eos->nz[i] == 2 && eos->na[i] == 4 ) xnuc[i] = xHe4;
    if(eos->nz[i] == 6 && eos->na[i] == 12) xnuc[i] = xC12;
    if(eos->nz[i] == 8 && eos->na[i] == 16) xnuc[i] = xO16;
    if(eos->nz[i] == 10 && eos->na[i] == 20) xnuc[i] = xNe20;
    if(eos->nz[i] == 12 && eos->na[i] == 24) xnuc[i] = xMg24;
    }

  eos_calc_tgiven( eos, rho0, xnuc, tempC, &res );
  rho[0] = rho0;
  r[0] = 1e2;
  dm[0] = rho0 * 4. / 3. * M_PI * r[0] * r[0] * r[0];
  p[0] = res.p.v;
  e[0] = res.e.v;
  temp[0] = tempC;
  He4[0] = xHeC;
  C12[0] = xCC;
  O16[0] = xOC;
  //printf("Rho: %e, p: %e, e: %e\n", rho[0], p[0], e[0]);

  double mass = dm[0];
  double rad = 1e2;
  double drad = 1e5;
  count = 1;
  var = 1;
  trans = 4; /* number of points in transition between core and shell */

  double y[2];
  y[0] = p[0];
  y[1] = 0.;

  double h[3];
  h[0] = y[0];
  h[1] = y[1];
  h[2] = temp[0];

  struct paramsWD params;
  params.temp = tempC;
  params.xnuc = xnuc;
  params.eos = eos;
  params.rho = rho0;


  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
        gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
        gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, 1e-6);
        gsl_odeiv_evolve * ev = gsl_odeiv_evolve_alloc (2);
        gsl_odeiv_system sys = {createWDIntegrator, 0, 2, &params};


  while (rho[count-1] > 1e-5)
    {
    if ((rho[count-1] >= rhoS)) /* correct for core */
      {
      int status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e10, &drad, y);

      if (status != GSL_SUCCESS)
        {
        PyErr_SetString( PyExc_ValueError, "ODE Solver failed.\n" );
        return 0;
        }

      /* increase arraylen if necessary */
      if (count == arraylen) {
        rho = resize( rho, arraylen, 2*arraylen );
        r = resize( r, arraylen, 2*arraylen );
        p = resize( p, arraylen, 2*arraylen );
        e = resize( e, arraylen, 2*arraylen );
        dm = resize( dm, arraylen, 2*arraylen );
        temp = resize( temp, arraylen, 2*arraylen );	
        He4 = resize( He4, arraylen, 2*arraylen );
        C12 = resize( C12, arraylen, 2*arraylen );
        O16 = resize( O16, arraylen, 2*arraylen );
        arraylen *= 2;
        }

      if(y[0] <= 0.)
        break;

      rho[count] = rho[count-1];

      eos_calc_ptgiven( eos, y[0], xnuc, tempC, &rho[count], &res );
      r[count] = rad;
      p[count] = y[0];
      e[count] = res.e.v;
      dm[count] = y[1] - mass;
      mass = y[1];
      temp[count] = temp[count-1];
      He4[count] = He4[count-1];
      C12[count] = C12[count-1];
      O16[count] = O16[count-1];

      }

    else if (rho[count-1] < rhoS) /* correct if HeShell */
      {

      if (var <= trans)
        {
        xHe4 = xHeS - (xHeS - xHeC) * (trans - var) / trans;
        xC12 = xCS - (xCS - xCC) * (trans - var) / trans;
        xO16 = xOS - (xOS - xOC) * (trans - var) / trans;
        He4[count] = xHe4;
        C12[count] = xC12;
        O16[count] = xO16;

        //temp[count] = tempS - 8 * deltaR * var; /* for adiab temp decrease in shell, but ony verz small decrease -> changed see below */
        //double xnuc[eos->nspecies];
        for(i = 0; i < eos->nspecies; i++)
          {
          xnuc[i] = 0.;
          if(eos->nz[i] == 2 && eos->na[i] == 4 ) xnuc[i] = xHe4;
          if(eos->nz[i] == 6 && eos->na[i] == 12) xnuc[i] = xC12;
          if(eos->nz[i] == 8 && eos->na[i] == 16) xnuc[i] = xO16;
          if(eos->nz[i] == 10 && eos->na[i] == 20) xnuc[i] = xNe20;
          if(eos->nz[i] == 12 && eos->na[i] == 24) xnuc[i] = xMg24;
          }

        params.xnuc = xnuc;

//        temp[count-1] = tempS;
        drad = 1e5;

        params.temp = temp[count-1];
        gsl_odeiv_system sys = {createWDIntegratorShell, 0, 3, &params};
        int status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e11, &drad, y);
        temp[count] = tempS - (tempS - tempC) * (trans - var) / trans;

        if (status != GSL_SUCCESS)
          {
          PyErr_SetString( PyExc_ValueError, "ODE Solver failed.\n" );
          return 0;
          }

          /* increase arraylen if necessary */
        if (count == arraylen) {
          rho = resize( rho, arraylen, 2*arraylen );
          r = resize( r, arraylen, 2*arraylen );
          p = resize( p, arraylen, 2*arraylen );
          e = resize( e, arraylen, 2*arraylen );
          dm = resize( dm, arraylen, 2*arraylen );
          temp = resize( temp, arraylen, 2*arraylen );
          He4 = resize( He4, arraylen, 2*arraylen );
          C12 = resize( C12, arraylen, 2*arraylen );
          O16 = resize( O16, arraylen, 2*arraylen );
          arraylen *= 2;
          }

        if(y[0] <= 0.){
          break;
          }
        rho[count] = rho[count-1];
        r[count] = rad;
        p[count] = y[0];

        eos_calc_ptgiven( eos, y[0], xnuc, temp[count], &rho[count], &res );
        e[count] = res.e.v;
        dm[count] = y[1] - mass;
        mass = y[1];
        }

      else
        {
        if (var == trans+1){
          h[0] = p[count-1];
          h[1] = mass;
          h[2] = temp[count-1];
          }

        xHe4 = xHeS;
        xC12 = xCS;
        xO16 = xOS;
        He4[count] = xHeS;
        C12[count] = xCS;
        O16[count] = xOS;

        //temp[count] = tempS - 8 * deltaR * var; /* for adiab temp decrease in shell, but ony verz small decrease -> changed see below */
        //double xnuc[eos->nspecies];
        for(i = 0; i < eos->nspecies; i++)
          {
          xnuc[i] = 0.;
          if(eos->nz[i] == 2 && eos->na[i] == 4 ) xnuc[i] = xHe4;
          if(eos->nz[i] == 6 && eos->na[i] == 12) xnuc[i] = xC12;
          if(eos->nz[i] == 8 && eos->na[i] == 16) xnuc[i] = xO16;
          if(eos->nz[i] == 10 && eos->na[i] == 20) xnuc[i] = xNe20;
          if(eos->nz[i] == 12 && eos->na[i] == 24) xnuc[i] = xMg24;
          }

        params.xnuc = xnuc;

        //params.temp = temp[count-1];
        gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 3);
        gsl_odeiv_evolve * ev = gsl_odeiv_evolve_alloc (3);
        gsl_odeiv_system sys = {createWDIntegratorShell, 0, 3, &params};
        int status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e11, &drad, h);
        drad = 1e6;

        if (status != GSL_SUCCESS)
          {
          PyErr_SetString( PyExc_ValueError, "ODE Solver failed.\n" );
          return 0;
          }

        /* increase arraylen if necessary */
        if (count == arraylen) {
          rho = resize( rho, arraylen, 2*arraylen );
          r = resize( r, arraylen, 2*arraylen );
          p = resize( p, arraylen, 2*arraylen );
          e = resize( e, arraylen, 2*arraylen );
          dm = resize( dm, arraylen, 2*arraylen );
          temp = resize( temp, arraylen, 2*arraylen );
          He4 = resize( He4, arraylen, 2*arraylen );
          C12 = resize( C12, arraylen, 2*arraylen );
          O16 = resize( O16, arraylen, 2*arraylen );
          arraylen *= 2;
          }

        if(h[0] <= 0.){
          break;
          }
        rho[count] = rho[count-1];
        r[count] = rad;
        p[count] = h[0];
        temp[count] = h[2];

        eos_calc_ptgiven( eos, h[0], xnuc, h[2], &rho[count], &res );
        e[count] = res.e.v;
        dm[count] = h[1] - mass;
        mass = h[1];
        }

      printf("%d, Radius: %e, Rho: %e, p: %e, e: %e, mass: %e\n", count, r[count], rho[count], p[count], e[count], mass);

      var++; /* for adiab temp decrease in shell */

      }

    if (rho[count] < 1e-5)
      break;
      
    if ((fabs(r[count]-r[count-1]) < 1e-5) && (fabs(r[count-2]-r[count-1]) < 1e-5)) /* to stop the run in case pressure and radius change only slightly due to bad initial conditions */
      break;

    if (h[2] <= 2e3){
      break;
      }
    count++;       

    }

  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, count ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, count ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, count ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, count ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, count ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( count ) );

  PyDict_SetStolenItem( dict, "temp", (PyObject*)createPyArray( temp, count ) );
  PyDict_SetStolenItem( dict, "He4", (PyObject*)createPyArray( He4, count ) );
  PyDict_SetStolenItem( dict, "C12", (PyObject*)createPyArray( C12, count ) );
  PyDict_SetStolenItem( dict, "O16", (PyObject*)createPyArray( O16, count ) );


  free( r );
  free( p );
  free( e );
  free( rho );
  free( dm );
  free( temp );
  free( He4 );
  free( C12 );
  free( O16 );

  gsl_odeiv_evolve_free(ev);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);


  return dict;
}

PyObject* _createWhiteDwarfHe(PyObject *self, PyObject *args) {
  double _rho, _temp; /* central density and temperature */
  double dr; /* width of radial zones */
  t_helm_eos_table *helm_eos_table;
  double xnuc[1];
  double *rho, *p, *e, *dm, *r;
  double rad, mass, dpdr;
  double ptmp, ptmp2, rho2;
  long arraylen;
  int i, n;
  int maxiter;
  PyObject* dict;
  struct eos_result res;

  maxiter = 20;
  dr = 1e4;
  if (!PyArg_ParseTuple( args, "O&dd|di:create_helmwdHe( eos, rho, temp, [dr, maxiter] )", &pyConvertHelmEos, &helm_eos_table, &_rho, &_temp, &dr, &maxiter )) {
    return 0;
  }
  
  arraylen = 0x10000;
  rho = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );
  
  xnuc[0] = 1.0;
  
  rho[0] = _rho;
  eos_calc_tgiven( helm_eos_table, rho[0], xnuc, _temp, &res );
  p[0] = res.p.v;
  e[0] = res.e.v;
  dm[0] = 0.0;
  
  mass = 0;
  rad = 0;
  i = 1;
  
  /* stop loop in case of integer overflow */
  while (i < 0x7FFFFFFF) {
    dpdr = -G*mass*rho[i-1] / max( rad*rad, 0.1 );
    rad += dr;
    p[i] = p[i-1] + dr*dpdr;

    rho[i] = rho[i-1];
    for (n=0; n<maxiter; n++) {
      eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
      ptmp = res.p.v;
      rho2 = 0.9999 * rho[i];
      eos_calc_tgiven( helm_eos_table, rho2, xnuc, _temp, &res );
      ptmp2 = res.p.v;
      rho[i] = rho[i] - (ptmp-p[i])/((ptmp-ptmp2)/(0.0001*rho[i]));
    }

    if (rho[i] < 1.0e-3) {
      printf( "EXIT: %g %g\n", rho[i], rad );
      break;
    }

    eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
    e[i] = res.e.v;
    ptmp = res.p.v;

    dm[i] = 0.5 * (rho[i-1]+rho[i]) * 4.0 * M_PI * pow( (rad-0.5*dr), 2.0 ) * dr;
    mass += dm[i];
    i++;
    

    if (i == arraylen) {
      rho = resize( rho, arraylen, 2*arraylen );
      p = resize( p, arraylen, 2*arraylen );
      e = resize( e, arraylen, 2*arraylen );
      dm = resize( dm, arraylen, 2*arraylen );
      arraylen *= 2;
    }
  }
  
  r = (double*)malloc( i * sizeof(double) );
  for (n=0; n<i; n++) {
    r[n] = dr*n;
  }
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, i ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, i ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, i ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, i ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, i ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( i ) );
  PyDict_SetStolenItem( dict, "dr", (PyObject*)PyFloat_FromDouble( dr ) );
  
  free( r );
  free( dm );
  free( p );
  free( e );
  free( rho );
 
  return dict;
}

PyObject* _createWhiteDwarfCO(PyObject *self, PyObject *args) {
  double _rho, _temp; /* central density and temperature */
  double xC12, xO16, xNe22; /* mass fractions */
  double dr; /* width of radial zones */
  t_helm_eos_table *helm_eos_table;
  double xnuc[100];
  double *rho, *p, *e, *dm, *r, *pot, *acc, *csnd;
  double rad, mass, dpdr;
  double ptmp, ptmp2, rho2;
  long arraylen;
  int i, n;
  int maxiter;
  PyObject* dict;
  struct eos_result res;

  maxiter = 20;
  dr = 1e3;
  xO16 = xNe22 = 0;
  if (!PyArg_ParseTuple( args, "O&ddd|dddi:create_wdCO( eos, rho, temp, xC12, [xO16, xNe22, dr, maxiter] )", &pyConvertHelmEos, &helm_eos_table, &_rho, &_temp, &xC12, &xO16, &xNe22, &dr, &maxiter )) {
    return 0;
  }
  
  arraylen = 0x10000;
  rho = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );
  pot = (double*)malloc( arraylen * sizeof(double) );
  acc = (double*)malloc( arraylen * sizeof(double) );
  csnd = (double*)malloc( arraylen * sizeof(double) );
  
  if (xO16 == 0. && xNe22 == 0.) {
    xnuc[0] = xC12;
    xnuc[1] = 1.0 - xO16;

    if (helm_eos_table->nspecies)
      xnuc[2] = 0.0;
  } else {
    xnuc[0] = xC12;
    xnuc[1] = xO16;
    xnuc[2] = xNe22;
    printf("Abundances: C12=%g O16=%g Ne22=%g, sum=%g.\n", xnuc[0], xnuc[1], xnuc[2], xnuc[0] + xnuc[1] + xnuc[2] );
  }
  
  rho[0] = _rho;
  eos_calc_tgiven( helm_eos_table, rho[0], xnuc, _temp, &res );
  e[0] = res.e.v;
  p[0] = res.p.v;
  dm[0] = 0.0;
  pot[0] = 0.0;
  csnd[0] = res.sound;
  
  mass = 0;
  rad = 0;
  i = 1;
  
  /* stop loop in case of integer overflow */
  while (i < 0x7FFFFFFF) {
    dpdr = -G*mass*rho[i-1] / max( rad*rad, 0.1 );
    acc[i] = dpdr;
    pot[i] = -G*mass / max(rad+0.5*dr, 0.1);
    rad += dr;
    p[i] = p[i-1] + dr*dpdr;

    rho[i] = rho[i-1];
    for (n=0; n<maxiter; n++) {
      eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
      ptmp = res.p.v;
      rho2 = 0.9999 * rho[i];
      eos_calc_tgiven( helm_eos_table, rho2, xnuc, _temp, &res );
      ptmp2 = res.p.v;
      rho[i] = rho[i] - (ptmp-p[i])/((ptmp-ptmp2)/(0.0001*rho[i]));
    }

    if (rho[i] < 1.0e-3) {
      printf( "EXIT: %g %g\n", rho[i], rad );
      break;
    }

    eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
    e[i] = res.e.v;
    ptmp = res.p.v;
    csnd[i] = res.sound;

    dm[i] = 0.5 * (rho[i-1]+rho[i]) * 4.0 * M_PI * pow( (rad-0.5*dr), 2.0 ) * dr;
    mass += dm[i];
    i++;
    

    if (i == arraylen) {
      rho = resize( rho, arraylen, 2*arraylen );
      p = resize( p, arraylen, 2*arraylen );
      e = resize( e, arraylen, 2*arraylen );
      dm = resize( dm, arraylen, 2*arraylen );
      pot = resize( pot, arraylen, 2*arraylen );
      acc = resize( acc, arraylen, 2*arraylen );
      csnd = resize( csnd, arraylen, 2*arraylen );
      arraylen *= 2;
    }
  }
  
  r = (double*)malloc( i * sizeof(double) );
  for (n=0; n<i; n++) {
    r[n] = dr*n;
  }
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, i ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, i ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, i ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, i ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, i ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( i ) );
  PyDict_SetStolenItem( dict, "dr", (PyObject*)PyFloat_FromDouble( dr ) );
  PyDict_SetStolenItem( dict, "pot", (PyObject*)createPyArray( pot, i ) );
  PyDict_SetStolenItem( dict, "acc", (PyObject*)createPyArray( acc, i ) );
  PyDict_SetStolenItem( dict, "csnd", (PyObject*)createPyArray( csnd, i ) );
  
  free( r );
  free( dm );
  free( p );
  free( e );
  free( rho );
  free( pot );
  free( acc );
  free( csnd );

  return dict;
}

PyObject* _createWhiteDwarfONeMg(PyObject *self, PyObject *args) {
  double _rho, _temp; /* central density and temperature */
  double _xO, _xNe, _xMg; /* mass fractions of oxygen, neon and magnesium */
  double dr; /* width of radial zones */
  t_helm_eos_table *helm_eos_table;
  double xnuc[3];
  double *rho, *p, *e, *dm, *r;
  double rad, mass, dpdr;
  double ptmp, ptmp2, rho2;
  long arraylen;
  int i, n;
  int maxiter;
  PyObject* dict;
  struct eos_result res;

  maxiter = 20;
  dr = 1e3;
  if (!PyArg_ParseTuple( args, "O&ddddd|di:create_wdONeMg( eos, rho, temp, xO, xNe, xMg, [dr, maxiter] )", &pyConvertHelmEos, &helm_eos_table, &_rho, &_temp, &_xO, &_xNe, &_xMg, &dr, &maxiter )) {
    return 0;
  }
  
  arraylen = 0x10000;
  rho = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );
  
  xnuc[0] = _xO;
  xnuc[1] = _xNe;
  xnuc[2] = _xMg;
  
  rho[0] = _rho;
  eos_calc_tgiven( helm_eos_table, rho[0], xnuc, _temp, &res );
  p[0] = res.p.v;
  e[0] = res.e.v;
  dm[0] = 0.0;
  
  mass = 0;
  rad = 0;
  i = 1;
  
  /* stop loop in case of integer overflow */
  while (i < 0x7FFFFFFF) {
    dpdr = -G*mass*rho[i-1] / max( rad*rad, 0.1 );
    rad += dr;
    p[i] = p[i-1] + dr*dpdr;

    rho[i] = rho[i-1];
    for (n=0; n<maxiter; n++) {
      eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
      ptmp = res.p.v;
      rho2 = 0.9999 * rho[i];
      eos_calc_tgiven( helm_eos_table, rho2, xnuc, _temp, &res );
      ptmp2 = res.p.v;
      rho[i] = rho[i] - (ptmp-p[i])/((ptmp-ptmp2)/(0.0001*rho[i]));
    }

    if (rho[i] < 1.0e-3) {
      printf( "EXIT: %g %g\n", rho[i], rad );
      break;
    }

    eos_calc_tgiven( helm_eos_table, rho[i], xnuc, _temp, &res );
    e[i] = res.e.v;
    ptmp = res.p.v;

    dm[i] = 0.5 * (rho[i-1]+rho[i]) * 4.0 * M_PI * pow( (rad-0.5*dr), 2.0 ) * dr;
    mass += dm[i];
    i++;
    

    if (i == arraylen) {
      rho = resize( rho, arraylen, 2*arraylen );
      p = resize( p, arraylen, 2*arraylen );
      e = resize( e, arraylen, 2*arraylen );
      dm = resize( dm, arraylen, 2*arraylen );
      arraylen *= 2;
    }
  }
  
  r = (double*)malloc( i * sizeof(double) );
  for (n=0; n<i; n++) {
    r[n] = dr*n;
  }
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, i ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, i ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, i ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, i ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, i ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( i ) );
  
  free( r );
  free( dm );
  free( p );
  free( e );
  free( rho );
  
  return dict;
}

PyObject* _createPolytrope(PyObject *self, PyObject *args, PyObject *keywds) {
  t_helm_eos_table *helm_eos_table;
  double rho0, pres0, temp0, K, dr;
  int maxiter;
  PyArrayObject *pyXnuc;
  double *xnuc, n, fac, ifac;
  double *rho, *p, *e, *dm, *r;
  double rad, mass, dpdr, temp;
  long arraylen;
  int i, j;
  PyObject* dict;
  struct eos_result res;

  static char *kwlist[] = { "eos", "n", "rho0", "composition", "pres0", "temp0", "dr", "maxiter", NULL };

  maxiter = 20;
  dr = 1e5;
  pres0 = 0.;
  temp0 = 0.;
  if (!PyArg_ParseTupleAndKeywords( args, keywds, "O&ddO!|dddi:create_polytrope( eos, n, rho0, composition, [pres0, temp0, dr, maxiter] )", kwlist, &pyConvertHelmEos, &helm_eos_table, &n, &rho0, &PyArray_Type, &pyXnuc, &pres0, &temp0, &dr, &maxiter )) {
    return 0;
  }

  xnuc = (double*)PyArray_DATA(pyXnuc);
  
  arraylen = 0x10000;
  rho = (double*)malloc( arraylen * sizeof(double) );
  p = (double*)malloc( arraylen * sizeof(double) );
  e = (double*)malloc( arraylen * sizeof(double) );
  dm = (double*)malloc( arraylen * sizeof(double) );

  fac = (n+1.) / n;
  ifac = 1. / fac;
  
  rho[0] = rho0;
  dm[0] = 0.0;

  temp = -1.;
  if (pres0 > 0) {
    p[0] = pres0;
    eos_calc_pgiven( helm_eos_table, rho[0], xnuc, p[0], &temp, &res );
    e[0] = res.e.v;
    K = p[0] / pow( rho[0], fac );
  } else if (temp0 > 0) {
    temp = temp0;
    eos_calc_tgiven( helm_eos_table, rho[0], xnuc, temp, &res );
    p[0] = res.p.v;
    e[0] = res.e.v;
    K = p[0] / pow( rho[0], fac );
  } else {
    PyErr_SetString( PyExc_ValueError, "pres0 or temp0 have to be set.\n" );
    return 0;
  }

  printf( "Central density: %g\n", rho[0] );
  printf( "Central temperature: %g\n", temp );
  printf( "Polytropic index: %g\n", n );
  printf( "Polytropic constant: %g\n", K );
  
  mass = 0;
  rad = 0;
  i = 1;
  
  /* stop loop in case of integer overflow */
  while (i < 0x7FFFFFFF) {
    dpdr = -G*mass*rho[i-1] / max( rad*rad, 0.1 );
    rad += dr;
    p[i] = p[i-1] + dr*dpdr;
    rho[i] = pow( p[i] / K, ifac );

    eos_calc_pgiven( helm_eos_table, rho[i], xnuc, p[i], &temp, &res );
    e[i] = res.e.v;

    dm[i] = 0.5 * (rho[i-1]+rho[i]) * 4.0 * M_PI * pow( (rad-0.5*dr), 2.0 ) * dr;
    mass += dm[i];

    if (rho[i] < 1.0e-3) {
      printf( "EXIT: %g %g\n", rho[i], rad );
      break;
    }

    i++;
    
    if (i == arraylen) {
      rho = resize( rho, arraylen, 2*arraylen );
      p = resize( p, arraylen, 2*arraylen );
      e = resize( e, arraylen, 2*arraylen );
      dm = resize( dm, arraylen, 2*arraylen );
      arraylen *= 2;
    }
  }
  
  r = (double*)malloc( i * sizeof(double) );
  for (j=0; j<i; j++) {
    r[j] = dr*j;
  }
  
  dict = PyDict_New();
  PyDict_SetStolenItem( dict, "rho", (PyObject*)createPyArray( rho, i ) );
  PyDict_SetStolenItem( dict, "u", (PyObject*)createPyArray( e, i ) );
  PyDict_SetStolenItem( dict, "p", (PyObject*)createPyArray( p, i ) );
  PyDict_SetStolenItem( dict, "dm", (PyObject*)createPyArray( dm, i ) );
  PyDict_SetStolenItem( dict, "r", (PyObject*)createPyArray( r, i ) );
  PyDict_SetStolenItem( dict, "ncells", (PyObject*)PyLong_FromLong( i ) );
  PyDict_SetStolenItem( dict, "dr", (PyObject*)PyFloat_FromDouble( dr ) );
  
  free( r );
  free( dm );
  free( p );
  free( e );
  free( rho );
  
  return dict;
}

static PyMethodDef icmethods[] = {
  { "create_wd", (PyCFunction)_createWhiteDwarf, METH_VARARGS|METH_KEYWORDS, "" },
  { "create_wdHeShell", (PyCFunction)_createWhiteDwarfHeShell, METH_VARARGS|METH_KEYWORDS, ""},
  { "create_wdHe", _createWhiteDwarfHe, METH_VARARGS, "" },
  { "create_wdCO", _createWhiteDwarfCO, METH_VARARGS, "" },
  { "create_wdONeMg", _createWhiteDwarfONeMg, METH_VARARGS, "" },
  { "create_polytrope", (PyCFunction)_createPolytrope, METH_VARARGS|METH_KEYWORDS, "" },
  { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "ic", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  icmethods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_ic(void)
{
	import_array();
	return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initic(void)
{
	Py_InitModule( "ic", icmethods );
	import_array();
}
#endif
