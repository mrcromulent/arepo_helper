#include <Python.h>
#include <arrayobject.h>
#include <string>
#include <cmath>
#include <random>

#include "helm_eos.h"
#include "sph.h"
#include "ic.h"
#include "main.h"
#include "const.h"
#include "utils.h"
#include "visualise.h"
#include "create_ics.h"
#include "write_ics.h"

int main() {

    if(PyArray_API == nullptr)
    {
        Py_Initialize();
        import_array()
    }

    test_eos();
    test_tree();
    test_tree_2();
    test_make_pcolor();
    test_make_radial();
    test_make_wd();
    test_make_polytrope();
    test_make_wdec_newest();
    test_make_he_wd();
    test_add_grid_particles();

    return 0;
}

void test_add_grid_particles() {
    const char *wdec_dir = "/home/pierre/wdec/";
    int nspecies = 5;
    double boxsize = 1e10;
    double pmass = 1e-6 * msol;
    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);


    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, 1.4);

    auto wd = convert_to_healpix_implementation(dict1, boxsize, centers, 0, 0, pmass);

    PyObject *xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);
    auto ret_dict = add_grid_particles_implementation(wd, boxsize, 32, 2e6, 1e-4, xnuc);
    print((PyObject *) ret_dict);

}

void test_make_he_wd() {
    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";
    double pmass                  = 1e-6 * msol;
    double mtot                   = 0.35; // msol
    double boxsize                = 1e10;
    const char* out_name          = "0_35He.hdf5";

    eos_init(eos, datafile, speciesfile, 0, 1);

    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);
    PyObject *xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);


    double temp_c = 5e5;
    double rho_c = rho_c_from_mtot_implementation(mtot, temp_c, eos, xnuc);

    auto dict = create_wd_implementation(eos, rho_c, temp_c, xnuc, 1e-6);
    eos_deinit(eos);

    auto dict2 = convert_to_healpix_implementation(dict, boxsize, centers, 0, 0, pmass);
    auto dict3 = add_grid_particles_implementation(dict2, boxsize, 32, 2e6, 1e-4, xnuc);

    print(PyDict_GetItemString(dict3, f[N::DENSITY]));
    print(PyDict_GetItemString(dict3, f[N::INTERNALENERGY]));

    write_dict_to_hdf5(dict3, out_name, boxsize);
}

void test_make_wdec_newest() {
    const char *wdec_dir = "/home/pierre/Desktop/WDEC/wdec0_75COHe/"; // TODO: Trailing slash is required!
    const char *out_name = "wdec0_75COHe.hdf5";
    double boxsize = 1e10;
    double pmass = 1e-6 * msol;
    int nspecies = 5;
    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);
    PyObject *xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);


    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, 1.4);

    // Print results
    print(PyDict_GetItemString(dict1, f[N::DENSITY]));
    print(PyDict_GetItemString(dict1, f[N::NUCLEARCOMPOSITION]));
    print(PyDict_GetItemString(dict1, f[N::RADIUS]));


    auto dict2 = convert_to_healpix_implementation(dict1, boxsize, centers, 0, 0, pmass);
    auto dict3 = add_grid_particles_implementation(dict2, boxsize, 32, 2e6, 1e-4, xnuc);

    print(PyDict_GetItemString(dict3, f[N::DENSITY]));
    print(PyDict_GetItemString(dict1, f[N::INTERNALENERGY]));


    write_dict_to_hdf5(dict3, out_name, boxsize);

}

void test_make_polytrope() {

    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 0, 1);


    PyObject *xnuc = Py_BuildValue("[ddddd]", 0.0, 0.5, 0.5, 0.0, 0.0);

    auto dict = create_polytrope_implementation(eos, 3, 5e6, xnuc, 0.0, 5e5, 1e6);
    eos_deinit(eos);

    print(PyDict_GetItemString(dict, f[N::DENSITY]));
    print(PyDict_GetItemString(dict, f[N::RADIUS]));
}

void test_make_wd() {
    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 0, 1);

    PyObject *xnuc = Py_BuildValue("[ddddd]", 0.0, 0.5, 0.5, 0.0, 0.0);

    auto dict = create_wd_implementation(eos, 5e6, 5e5, xnuc, 1e-6);
    eos_deinit(eos);

}

void test_make_pcolor() {

    const char *wdec_dir = "/home/pierre/wdec/";
    int nspecies = 5;
    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, 1.4);

    double boxsize = 1e10;
    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);



    auto grid_xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);
    auto wd = convert_to_healpix_implementation(dict1, boxsize, centers,0, 0, 1e-6 * msol);
    auto wd2 = add_grid_particles_implementation(wd, boxsize, 32, 2e6, 1e-4, grid_xnuc);

    auto value = (PyArrayObject *) PyDict_GetItemString(wd2, f[N::DENSITY]);
    auto pos = (PyArrayObject *) PyDict_GetItemString(wd2, f[N::COORDINATES]);
    print((PyObject *) value);
    print((PyObject *) pos);


    PyObject *boxsizes      = Py_BuildValue("[dd]", boxsize, boxsize);
    PyObject *resolutions   = Py_BuildValue("[ii]", 1024, 1024);
    PyObject *axes          = Py_BuildValue("[ii]", 0, 1);
    int include_neighbours_in_output = 1;
    int numthreads          = 1;

    auto dict = make_pcolor_implementation(pos, value, axes, boxsizes, resolutions, centers,
                                           include_neighbours_in_output,
                                           numthreads);

    print(PyDict_GetItemString(dict, "grid"));
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
    auto *xnuc = (double *) malloc(nspecies * sizeof(double));
    for (int i = 0; i < nspecies; i++) {
        xnuc[i] = 0;
    }

    xnuc[0] = 0.5;
    xnuc[1] = 0.5;
    double tempguess = 1e8;
    eos_result res{};
    eos_calc_pgiven(helm_eos_table, 1e6, xnuc, 1e22, &tempguess, &res);

    for (int i = 0; i < nspecies; i++) {
        xnuc[i] = 0;
    }
    xnuc[0] = 1.0;
    double tempguess2 = 1500;
    eos_result res2{};
    eos_calc_egiven(helm_eos_table, 1e-7, xnuc, 1e13, &tempguess2, &res2);
    std::cout << res.p.v << std::endl;

    eos_deinit(helm_eos_table);

}

void test_tree() {
    t_sph_tree tree;

    int npart = 2;
    auto *pos = (double *) malloc(3 * npart * sizeof(double));

    pos[0 * 3 + 0] = 0;
    pos[0 * 3 + 1] = 0;
    pos[0 * 3 + 2] = 0;
    pos[1 * 3 + 0] = 1;
    pos[1 * 3 + 1] = 1;
    pos[1 * 3 + 2] = 1;

    createTree(&tree, npart, pos);
    // -------------------------------------------------------------------------------------------------------------
}

void test_tree_2() {

    const char *wdec_dir = "/home/pierre/wdec/";
    int nspecies = 5;
    double boxsize = 1e10;
    double gamma = 1.4;
    double pmass = 1e-6 * msol;
    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);

    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, gamma);

    auto wd = convert_to_healpix_implementation(dict1, boxsize,
                                                centers, 0, 0, pmass);

    auto density_1 = (PyArrayObject *) PyDict_GetItemString(wd, f[N::DENSITY]);
    auto pos_1 = (PyArrayObject *) PyDict_GetItemString(wd, f[N::COORDINATES]);
    auto mass_1 = (PyArrayObject *) PyDict_GetItemString(wd, f[N::MASSES]);
    int npart = PyArray_DIMS(density_1)[0];
    auto pos = (double *) PyArray_DATA(pos_1);
    auto density = (double *) PyArray_DATA(density_1);
    auto mass = (double *) PyArray_DATA(mass_1);

    t_sph_tree tree;
    createTree(&tree, npart, pos);

    // Coordinate
    auto *coord = (double *) malloc(3 * sizeof(double));
    coord[0] = boxsize/2; coord[1] = boxsize/2; coord[2] = boxsize/2;

    // Domain length
    auto domainLen = getDomainLen(npart, pos);

    // Get Nearest Node/Particles
    auto node = getNearestNode(&tree, coord);
    auto particles = getParticles(&tree, node);

    double hsml=0;
    double *dhsmldensity;
    auto weighted_neighbours = calcHsml(&tree, coord, pos, mass, 4, &hsml, density);
    auto weighted_neighbours_2 = calcDensity(&tree, coord, hsml, pos, mass, density,
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

void test_make_radial() {

    const char *wdec_dir    = "/home/pierre/wdec/";
    int nspecies            = 5;
    double pmass            = 1e-6 * msol;
    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, 1.4);

    double boxsize      = 1e10;
    PyObject *centers   = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);
    PyObject *b         = Py_BuildValue("[ddd]", boxsize, boxsize/2, boxsize/2);
    auto a              = centers;
    double cylinder_rad = 0.01 * boxsize;
    int nshells         = 200;

    auto wd = convert_to_healpix_implementation(dict1, boxsize, centers, 0, 0, pmass);

    auto pos = (PyArrayObject *) PyDict_GetItemString(wd, f[N::COORDINATES]);
    auto value = (PyArrayObject *) PyDict_GetItemString(wd, f[N::DENSITY]);


    auto radial = make_radial_implementation(pos, value, a, b, cylinder_rad, nshells);
    print((PyObject *) radial);

}
