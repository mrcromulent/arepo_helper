#include <Python.h>
#include <hdf5.h>
#include <arrayobject.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <cmath>

#include <h5cpp/core>
#include <h5cpp/io>

#include "helm_eos.h"
#include "sph.h"
#include "ic.h"
#include "main.h"
#include "omp_util.h"
#include "const.h"
#include "make_wdec.h"
#include "utils.h"
#include "make_pcolor.h"

void test_make_pcolor();

int main() {

    if(PyArray_API == nullptr)
    {
        Py_Initialize();
        import_array()
    }


//    test_eos();
//    test_tree();
//    test_WD();
//    test_tree_2();
//    test_make_wdec();
    test_make_pcolor();

    return 0;
}

void test_make_pcolor() {

    const char *wdec_dir    = "/home/pierre/wdec/";
    double pmass            = 1e-6 * msol;

    auto wd = make_wdec_new(wdec_dir, 1e10, 5, true, false, false,
                            pmass);
    auto value = (PyArrayObject *) PyDict_GetItemString(wd, "rho");
    auto pos = (PyArrayObject *) PyDict_GetItemString(wd, "pos");

//    auto value_data = (double *) PyArray_DATA(value);
//    auto value_pos = (double *) PyArray_DATA(pos);
//    auto n_value = PyArray_DIMS(value)[0];
//
//    std::cout << n_value << std::endl;
//    for (int i = 0; i < n_value; i++) {
//        if (i % 10000 == 0) {
//            std::cout << " i = " << i << ". x = " << value_pos[i * 3 + 0] << ". q = " << value_data[i] << std::endl;
//        }
//    }

    double boxsize = 1e10;
    int numthreads = 1;
    double boxsize_x = boxsize, boxsize_y = boxsize, boxsize_z = 0;
    double center_x = boxsize/2, center_y = boxsize/2, center_z = boxsize/2;
    int resolution_x = 1024;
    int resolution_y = 1024;
    // These dictate what is on the x and y axis
    int axis0 = 0;
    int axis1 = 1;
    int nz = 1;
    bool include_neighbours_in_output = true;

    auto dict = make_pcolor_implementation(pos, value,
                                           resolution_x, resolution_y,
                                           boxsize_x, boxsize_y, boxsize_z,
                                           center_x, center_y, center_z,
                                           axis0, axis1,
                                           nz,
                                           include_neighbours_in_output,
                                           numthreads);

    auto grid = (PyArrayObject *) PyDict_GetItemString(dict, "grid");
    auto grid_data = (double *) PyArray_DATA(grid);
    auto n = PyArray_DIMS(grid)[0];

    std::cout << n << std::endl;
    for (int i = 0; i < n; i++) {
        if (i % 10 == 0) {
            std::cout << grid_data[i] << std::endl;
        }
    }
    std::cout << std::endl;
}


void test_eos() {
    // LOADING A HELM EOS FILE
    // -------------------------------------------------------------------------------------------------------------
    auto *helm_eos_table = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(helm_eos_table, datafile, speciesfile, 0, 1);

    // Calling the 'given' functions
    // -------------------------------------------------------------------------------------------------------------

    int nspecies = 5;
    auto *xnuc = static_cast<double *>(malloc(nspecies * sizeof(double)));

    { int i;
        for (i=0; i<nspecies; i++) xnuc[i] = 0;
    }

    xnuc[0] = 0.5;
    xnuc[1] = 0.5;
    double tempguess = 1e8;
    eos_result res{};
    eos_calc_pgiven(helm_eos_table, 1e6, xnuc, 1e22, &tempguess, &res);

    eos_deinit(helm_eos_table);
    // -------------------------------------------------------------------------------------------------------------
}

void test_tree() {
    // CREATE A TREE
    // -------------------------------------------------------------------------------------------------------------
    t_sph_tree tree;

    int npart = 2;
    auto *pos = (double *) malloc(3 * npart * sizeof(double));

    pos[0 * 3 + 0] = 0;
    pos[0 * 3 + 1] = 0;
    pos[0 * 3 + 2] = 0;
    pos[1 * 3 + 0] = 1;
    pos[1 * 3 + 1] = 1;
    pos[1 * 3 + 2] = 1;

    createTree( &tree, npart, pos );
    // -------------------------------------------------------------------------------------------------------------
}

void test_tree_2() {

    const char *wdec_dir    = "/home/pierre/wdec/";
    double pmass            = 1e-6 * msol;

    auto wd = make_wdec_new(wdec_dir, 1e10, 5, true, false, false,
                            pmass);

    auto density_1 = (PyArrayObject *) PyDict_GetItemString(wd, "rho");
    auto pos_1 = (PyArrayObject *) PyDict_GetItemString(wd, "pos");
    auto mass_1 = (PyArrayObject *) PyDict_GetItemString(wd, "mass");
    int npart = PyArray_DIMS(density_1)[0];
    auto pos = (double *) PyArray_DATA(pos_1);
    auto density = (double *) PyArray_DATA(density_1);
    auto mass = (double *) PyArray_DATA(mass_1);

    t_sph_tree tree;
    createTree( &tree, npart, pos );

    // Coordinate
    auto *coord = (double *) malloc(3 * sizeof(double));
    coord[0] = 1e10/2; coord[1] = 1e10/2; coord[2] = 1e10/2;

    // Domain length
    auto domainLen = getDomainLen(npart, pos);

    // Get Nearest Node/Particles
    auto node = getNearestNode(&tree, coord);
    auto particles = getParticles(&tree, node);

    double hsml=0;
    double weighted_neighbours, weighted_neighbours_2;
    double *dhsmldensity;
    weighted_neighbours = calcHsml(&tree, coord, pos, mass, 4, &hsml, density);
    weighted_neighbours_2 = calcDensity(&tree, coord, hsml, pos, mass, density,
                                        reinterpret_cast<double *>(&dhsmldensity));

    //
    int nneighbours_real = 4;
    int *neighbours = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));
    int converged;
    auto nneighbours = getNeighbours(&tree, coord, pos, hsml, reinterpret_cast<int **>(neighbours));

    auto radius = getNNeighbours(&tree, coord, pos, 4,
                                 &nneighbours_real,
                                 reinterpret_cast<int **>(&neighbours),
                                 &converged);

    // Free tree
    freeTree(&tree);
}

void test_make_wdec() {
    const char *wdec_dir    = "/home/pierre/wdec/";
    double pmass            = 1e-6 * msol;

    auto wd = make_wdec_new(wdec_dir, 1e10, 5, true, false, false,
                            pmass);
    auto value = (PyArrayObject *) PyDict_GetItemString(wd, "rho");
    auto pos = (PyArrayObject *) PyDict_GetItemString(wd, "pos");

    auto value_data = (double *) PyArray_DATA(value);
    auto value_pos = (double *) PyArray_DATA(pos);
    auto n_value = PyArray_DIMS(value)[0];

    std::cout << n_value << std::endl;
    for (int i = 0; i < n_value; i++) {
        if (i % 10000 == 0) {
            std::cout << " i = " << i << ". x = " << value_pos[i * 3 + 0] << ". q = " << value_data[i] << std::endl;
        }
    }
}

PyObject* test_WD() {
    double rho0, temp; /* central density and temperature */
    double xHe4, xC12, xO16, xNe20, xMg24; /* mass fractions */
    double tol;
//    t_helm_eos_table *eos;
    double *r, *p, *e, *rho, *dm, *csnd;
    int i, arraylen, count;
    PyObject* dict;
    struct eos_result res = {};

    // Set these
    auto *eos   = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 0, 0);
    tol = 1e-6;
    rho0 = 1e7;
    temp = 5e5;
    xHe4 = 0.0;
    xC12 = 0.5;
    xO16 = 0.5;
    xNe20 = xMg24 = 0.;

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

    struct paramsWD params = {};
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
    ///
    PyDict_SetStolenItem( dict, "xnuc", (PyObject*)createPyArray( xnuc,  eos->nspecies) );
    ///
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
