#include <Python.h>
#include <arrayobject.h>
#include <cmath>

#include "headers/helm_eos.h"
#include "headers/sph.h"
#include "headers/ic.h"
#include "main.h"
#include "headers/const.h"
#include "headers/utils.h"
#include "headers/visualise.h"
#include "headers/create_ics.h"
#include "headers/write_ics.h"

int main() {

    if(PyArray_API == nullptr)
    {
        Py_Initialize();
        import_array()
    }

//    test_eos();
//    test_tree();
//    test_tree_2();
//    test_make_pcolor();
//    test_make_radial();
    test_make_wd();
//    test_make_polytrope();
//    test_make_wdec_newest();
//    test_make_he_wd();
//    test_add_grid_particles();

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
    auto ret_dict = add_grid_particles_implementation(wd, boxsize, 32, 4e6, 1e-4, xnuc);
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

    eos_init(eos, datafile, speciesfile, 1, 1);

    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);
    PyObject *xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);


    double temp_c = 5e5;
    double rho_c = rho_c_from_mtot_implementation(mtot, temp_c, eos, xnuc);

    auto dict = create_wd_implementation(eos, rho_c, temp_c, xnuc, 1e-6);
    eos_deinit(eos);

    auto dict2 = convert_to_healpix_implementation(dict, boxsize, centers, 0, 0, pmass);
    auto dict3 = add_grid_particles_implementation(dict2, boxsize, 32, 4e6, 1e-4, xnuc);

    print(dict3);

//    write_dict_to_hdf5(dict3, out_name, boxsize);
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
    auto dict3 = add_grid_particles_implementation(dict2, boxsize, 32, 4e6, 1e-4, xnuc);

    print(PyDict_GetItemString(dict3, f[N::DENSITY]));
    print(PyDict_GetItemString(dict1, f[N::INTERNALENERGY]));


    write_dict_to_hdf5(dict3, out_name, boxsize);

}

void test_make_polytrope() {

    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 1, 1);


    PyObject *xnuc = Py_BuildValue("[ddddd]", 0.0, 0.5, 0.5, 0.0, 0.0);

    auto dict = create_polytrope_implementation(eos, 3, 5e6, xnuc, 0.0, 5e5, 1e6);
    eos_deinit(eos);

    print(PyDict_GetItemString(dict, f[N::DENSITY]));
    print(PyDict_GetItemString(dict, f[N::RADIUS]));
}

void test_make_wd() {
    auto *eos = (t_helm_eos_table *) malloc(sizeof(t_helm_eos_table));
    const char *datafile          = "/home/pierre/Downloads/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(eos, datafile, speciesfile, 1, 0);

    PyObject *xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);

    double temp_c = 5e5;
    double rho_c = rho_c_from_mtot_implementation(0.35, temp_c, eos, xnuc);
    auto dict = create_wd_implementation(eos, rho_c, temp_c, xnuc, 1e-6);
//    print(dict);
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
    auto wd2 = add_grid_particles_implementation(wd, boxsize, 32, 4e6, 1e-4, grid_xnuc);

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
    const char *datafile          = "/home/pierre/Downloads/helm_table.dat";
    const char *speciesfile       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt";

    eos_init(helm_eos_table, datafile, speciesfile, 1, 1);

    // Calling the 'given' functions
    // -------------------------------------------------------------------------------------------------------------

    int nspecies = 5;
    auto *xnuc = (double *) malloc(nspecies * sizeof(double));
    for (int i = 0; i < nspecies; i++) {
        xnuc[i] = 0;
    }

    xnuc[0] = 0.5;
    xnuc[1] = 0.5;
    double tempguess = -1;
    eos_result res{};
    eos_calc_pgiven(helm_eos_table, 1e6, xnuc, 1e22, &tempguess, &res);

    for (int i = 0; i < nspecies; i++) {
        xnuc[i] = 0;
    }
    xnuc[0] = 1.0;
    double tempguess2 = -1;
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
    auto grid_xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);
    PyObject *centers = Py_BuildValue("[ddd]", boxsize/2, boxsize/2, boxsize/2);

    auto dict1 = create_wd_wdec_implementation(wdec_dir, nspecies, gamma);

    auto wd = convert_to_healpix_implementation(dict1, boxsize, centers, 0, 0, pmass);
    auto wd2 = add_grid_particles_implementation(wd, boxsize, 32, 4e6, 1e-4, grid_xnuc);

    auto density_1 = (PyArrayObject *) PyDict_GetItemString(wd2, f[N::DENSITY]);
    auto density = (double *) PyArray_DATA(density_1);
    int npart = PyArray_DIMS(density_1)[0];
    auto pos = (double *) PyArray_DATA((PyArrayObject *) PyDict_GetItemString(wd2, f[N::COORDINATES]));
    auto mass = (double *) PyArray_DATA((PyArrayObject *) PyDict_GetItemString(wd2, f[N::MASSES]));

    t_sph_tree tree;
    createTree(&tree, npart, pos);

    // Forward definitions
    double *coord;
    double hsml, x, y, z, domainLen, weighted_neighbours;
    int n_neighbours, n, node;
    int *neighbours;
    t_sph_treenode pnode;

    // Coordinate
    coord = (double *) malloc(3 * sizeof(double));

    // Number of neighbours is n_neighbours
    // Neighbour numbers stored in neighbours
    hsml = 1e7;
    neighbours = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));
    coord[0] = boxsize; coord[1] = boxsize; coord[2] = boxsize;
    n_neighbours = getNeighbours(&tree, coord, pos, hsml, &neighbours);

    for (int i = 0; i < n_neighbours; i++) {
        n = neighbours[i];
        if (n < tree.npart) { // particle node, that contains only one particle
            x = pos[n * 3 + 0];
            y = pos[n * 3 + 1];
            z = pos[n * 3 + 2];

        } else { // real node
            pnode = tree.nodes[n];
            x = pnode.center[0];
            y = pnode.center[1];
            z = pnode.center[2];
        }
        printf("n = %i; x, y, z = %g, %g, %g \n", n + 1, x, y, z);
    }

    // Domain length
    domainLen = getDomainLen(npart, pos);
    printf("domainLen = %g \n", domainLen);

    coord[0] = boxsize; coord[1] = boxsize; coord[2] = boxsize;
    node = getNearestNode(&tree, coord);
    if (node < tree.npart) { // particle node, that contains only one particle
        x = pos[node * 3 + 0];
        y = pos[node * 3 + 1];
        z = pos[node * 3 + 2];

    } else { // we found a real node
        pnode = tree.nodes[node];
        x = pnode.center[0];
        y = pnode.center[1];
        z = pnode.center[2];
    }
    printf("node = %i; x, y, z = %g, %g, %g \n", node + 1, x, y, z);

    // Calculates the smoothing length
    hsml = 0;
    n_neighbours = 4;
    coord[0] = boxsize / 2; coord[1] = boxsize / 2; coord[2] = boxsize / 2;
    weighted_neighbours = calcHsml(&tree, coord, pos, mass, n_neighbours, &hsml, density);
    printf("weighted_neighbours = %f. Nneighbours = %i \n", weighted_neighbours, n_neighbours);
    printf("Smoothing length = %g \n", hsml);


    double dhsmldensity = 0;
    hsml = 1e7;
    coord[0] = boxsize / 2; coord[1] = boxsize / 2; coord[2] = boxsize / 2;
    weighted_neighbours = calcDensity(&tree, coord, hsml, pos, mass, density, &dhsmldensity);
    printf("weighted_neighbours = %g. \n", weighted_neighbours);

    //
    n_neighbours = 4;
    neighbours = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));
    coord[0] = 0.2 * boxsize; coord[1] = 0.2 * boxsize; coord[2] = 0.2 * boxsize;
    int nneighbours_real;
    int converged;

    auto radius = getNNeighbours(&tree, coord, pos,
                                 n_neighbours,
                                 &nneighbours_real,
                                 &neighbours,
                                 &converged);

    printf("Radius = %g. converged = %i \n", radius, converged);

    for (int i = 0; i < nneighbours_real; i++) {
        n = neighbours[i];
        if (n < tree.npart) { // particle node, that contains only one particle
            x = pos[n * 3 + 0];
            y = pos[n * 3 + 1];
            z = pos[n * 3 + 2];

        } else { // real node
            pnode = tree.nodes[n];
            x = pnode.center[0];
            y = pnode.center[1];
            z = pnode.center[2];
        }
        printf("n = %i; x, y, z = %g, %g, %g \n", n + 1, x, y, z);
    }

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
