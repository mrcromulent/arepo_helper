#include <Python.h>
#include <arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "helm_eos.h"
#include "pyhelm_eos.h"
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "ic.h"
#include "const.h"
#include "utils.h"

int createWDIntegrator( double r, const double y[], double f[], void *params)
{
    double temp = ((struct paramsWD*)params)->temp;
    double *xnuc = ((struct paramsWD*)params)->xnuc;
    t_helm_eos_table *eos = ((struct paramsWD*)params)->eos;
    double rho = ((struct paramsWD*)params)->rho;

    struct eos_result res{};
    eos_calc_ptgiven( eos, y[0], xnuc, temp, &rho, &res );
    r = std::max( r, 1e2 );
    double r2 = r*r;

    f[0] = - G * y[1] * rho / r2;
    f[1] = 4. * M_PI  * r2 * rho;

    ((struct paramsWD*)params)->rho = rho;

    return GSL_SUCCESS;
}

PyObject* _createWhiteDwarf(PyObject *self, PyObject *args, PyObject *kwargs) {
    double rho0, temp; /* central density and temperature */
    double xHe4, xC12, xO16, xNe20, xMg24; /* mass fractions */
    double tol;
    t_helm_eos_table *eos;
    double *r, *p, *e, *rho, *dm, *csnd;
    int i, arraylen, count;
    PyObject* dict;
    struct eos_result res{};

    const char *kwlist[] = { "eos", "rho0", "temp", "xHe4", "xC12", "xO16", "xNe20", "xMg24", "tol", nullptr };
    auto keywords = (char **) kwlist;

    tol = 1e-6;
    temp = 5e5;
    xHe4 = xC12 = xO16 = xNe20 = xMg24 = 0.;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&d|ddddddd:create_wd( eos, rho0, [temp, xHe4, xC12, xO16, xNe20, xMg24, tol] )", keywords, &pyConvertHelmEos, &eos, &rho0, &temp,
                                     &xHe4, &xC12, &xO16, &xNe20, &xMg24, &tol )) {
        return nullptr;
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
        return nullptr;
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

    struct paramsWD params {};
    params.temp = temp;
    params.xnuc = xnuc;
    params.eos = eos;
    params.rho = rho0;

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.0, tol);
    gsl_odeiv_evolve * ev = gsl_odeiv_evolve_alloc (2);
    gsl_odeiv_system sys = {createWDIntegrator, nullptr, 2, &params};

    double mass = dm[0];
    double rad = 1e2;
    double drad = 1e3;
    while (rho[count-1] > 1e-5 && rad < 1e10)
    {
        int status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e10, &drad, y);

        if (status != GSL_SUCCESS)
        {
            PyErr_SetString( PyExc_ValueError, "ODE Solver failed.\n" );
            return nullptr;
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

PyObject* _createPolytrope(PyObject *self, PyObject *args, PyObject *kwargs) {
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
    struct eos_result res{};

    const char *kwlist[] = { "eos", "n", "rho0", "composition", "pres0", "temp0", "dr", "maxiter", nullptr };
    auto keywords = (char **) kwlist;

    maxiter = 20;
    dr = 1e5;
    pres0 = 0.;
    temp0 = 0.;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&ddO!|dddi:create_polytrope( eos, n, rho0, composition, [pres0, temp0, dr, maxiter] )", keywords, &pyConvertHelmEos, &helm_eos_table, &n, &rho0, &PyArray_Type, &pyXnuc, &pres0, &temp0, &dr, &maxiter )) {
        return nullptr;
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
        return nullptr;
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
        dpdr = -G*mass*rho[i-1] / std::max( rad*rad, 0.1 );
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
            auto curr_length    = (int) arraylen;
            auto new_length     = static_cast<int>(2*arraylen);

            rho = resize(rho, curr_length, new_length);
            p = resize(p, curr_length, new_length);
            e = resize(e, curr_length, new_length);
            dm = resize(dm, curr_length, new_length);
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

// Python Module Definition
static PyMethodDef icmethods[] = {
        { "create_wd", (PyCFunction)_createWhiteDwarf, METH_VARARGS|METH_KEYWORDS, "" },
        { "create_polytrope", (PyCFunction)_createPolytrope, METH_VARARGS|METH_KEYWORDS, "" },
        { nullptr, nullptr, 0, nullptr }
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "ic",
        nullptr,
        -1,
        icmethods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_ic(void)
{
    import_array();

    return PyModule_Create(&moduledef);
}
