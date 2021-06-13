#include <Python.h>
#include <arrayobject.h>
#include <cmath>

#include "helm_eos.h"
#include "sph.h"
#include "ic.h"
#include "main.h"
#include "omp_util.h"
#include "const.h"
#include "utils.h"
#include "make_pcolor.h"
#include "make_radial.h"
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


    return 0;
}

void test_make_wdec_newest() {
    const char *wdec_dir = "/home/pierre/wdec/";
    auto dict1 = create_wd_wdec_implementation(wdec_dir);
//    auto rho = (PyArrayObject *) PyDict_GetItemString(dict1, "rho");
//    auto r = (PyArrayObject *) PyDict_GetItemString(dict1, "r");
//    print_pyobj(r);
//    print_pyobj(rho);


    double boxsize = 1e10;
    double centers[3] = {boxsize/2, boxsize/2, boxsize/2};
    auto centers_py = (PyArrayObject *) createPyArray(centers, 3);

    auto dict2 = convert_to_healpix_implementation(dict1, 5, boxsize,
                                                centers_py, 1, 0, 0, 1e-6 * msol);

//    auto ndata_rho = (PyArrayObject *) PyDict_GetItemString(dict2, "ndata_rho");
//    auto ndata_u = (PyArrayObject *) PyDict_GetItemString(dict2, "ndata_u");
//    print_pyobj(ndata_rho);
//    print_pyobj(ndata_u);

    write_healpix_to_file(dict2, "bin.dat.ic.hdf5");

}

void test_make_polytrope() {

    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 0, 1);

    double xnuc[eos->nspecies];
    xnuc[0] = xnuc[3] = xnuc[4] = 0.0;
    xnuc[1] = xnuc[2] = 0.5;
    auto xnuc_py = (PyArrayObject *) createPyArray(xnuc, eos->nspecies);

    auto dict = create_polytrope_implementation(eos, 3, 5e6, xnuc_py, 0.0, 5e5, 1e6);
    eos_deinit(eos);

    auto rho = (PyArrayObject *) PyDict_GetItemString(dict, "rho");
    auto r = (PyArrayObject *) PyDict_GetItemString(dict, "r");
    print_pyobj(r);
    print_pyobj(rho);
}

void test_make_wd() {
    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 0, 1);

    double xnuc[eos->nspecies];
    xnuc[0] = xnuc[3] = xnuc[4] = 0.0;
    xnuc[1] = xnuc[2] = 0.5;
    auto xnuc_py = (PyArrayObject *) createPyArray(xnuc, eos->nspecies);

    auto dict = create_wd_implementation(eos, 5e6, 5e5, xnuc_py, 1e-6);
    eos_deinit(eos);

    auto rho = (PyArrayObject *) PyDict_GetItemString(dict, "rho");
    auto r = (PyArrayObject *) PyDict_GetItemString(dict, "r");
    print_pyobj(r);
    print_pyobj(rho);
}

void test_make_pcolor() {

    const char *wdec_dir = "/home/pierre/wdec/";
    auto dict1 = create_wd_wdec_implementation(wdec_dir);

    double boxsize = 1e10;
    double centers[3] = {boxsize/2, boxsize/2, boxsize/2};
    auto centers_py = (PyArrayObject *) createPyArray(centers, 3);


    auto wd = convert_to_healpix_implementation(dict1, 5, boxsize,
                                                centers_py, 1, 0, 0, 1e-6 * msol);

    auto value = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_rho");
    auto pos = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_pos");

    auto value_data = (double *) PyArray_DATA(value);
    auto value_pos = (double *) PyArray_DATA(pos);
    auto n_value = PyArray_DIMS(value)[0];

    std::cout << n_value << std::endl;
    for (int i = 0; i < n_value; i++) {
        if (i % 10000 == 0) {
            std::cout << " i = " << i << ". x = " << value_pos[i * 3 + 0] << ". q = " << value_data[i] << std::endl;
        }
    }

    double bs[2]        = {boxsize, boxsize};
    auto boxsizes       = (PyArrayObject *) createPyArray(bs, 2);
    int rs[2]           = {1024, 1025};
    auto resolutions    = (PyArrayObject *) createPyArray(rs, 2);
    int ax[2]           = {0, 1};
    auto axes           = (PyArrayObject *) createPyArray(ax, 2);




    int include_neighbours_in_output = 1;
    int numthreads = 1;

    auto dict = make_pcolor_implementation(pos, value, axes, boxsizes, resolutions, centers_py,
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

    const char *wdec_dir = "/home/pierre/wdec/";
    auto dict1 = create_wd_wdec_implementation(wdec_dir);

    double boxsize = 1e10;
    double centers[3] = {boxsize/2, boxsize/2, boxsize/2};
    auto centers_py = (PyArrayObject *) createPyArray(centers, 3);

    auto wd = convert_to_healpix_implementation(dict1, 5, boxsize,
                                                centers_py, 1, 0, 0, 1e-6 * msol);

    auto density_1 = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_rho");
    auto pos_1 = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_pos");
    auto mass_1 = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_mass");
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

void test_make_radial() {

    const char *wdec_dir = "/home/pierre/wdec/";
    auto dict1 = create_wd_wdec_implementation(wdec_dir);

    double boxsize      = 1e10;
    double centers[3]   = {boxsize/2, boxsize/2, boxsize/2};
    double b_val[3]     = {boxsize, boxsize/2, boxsize/2};
    auto centers_py     = (PyArrayObject *) createPyArray(centers, 3);
    auto b              = (PyArrayObject *) createPyArray(b_val, 3);
    auto a              = centers_py;
    double cylinder_rad = 0.01 * boxsize;

    auto wd = convert_to_healpix_implementation(dict1, 5, boxsize,
                                                centers_py, 1, 0, 0, 1e-6 * msol);

    auto pos = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_pos");
    auto value = (PyArrayObject *) PyDict_GetItemString(wd, "ndata_rho");
    int nshells = 200;

    auto returned_thing = make_radial_implementation(pos, value, a, b, cylinder_rad, nshells);

    auto value_data = (double *) PyArray_DATA(returned_thing);
    for (int i = 0; i < nshells; i++) {
        std::cout << value_data[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < nshells; i++) {
        std::cout << value_data[nshells + i] << " ";
    }
    std::cout << std::endl;

}
