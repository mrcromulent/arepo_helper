#include <Python.h>
#include <armadillo>
#include <arrayobject.h>

#include "../headers/const.h"
#include "../headers/sph.h"
#include "../headers/utils.h"
#include "../headers/omp_util.h"
#include "../headers/visualise.h"

// This uses the method described here:
// https://math.stackexchange.com/questions/3518495/check-if-a-general-point-is-inside-a-given-cylinder
// Note that this method fails if the point A or B is the origin!
bool inside_cylinder(const arma::vec& A, const arma::vec& B, const arma::vec& P, double radius, arma::vec& Q) {
    using namespace arma;

    vec E = B - A;
    vec m = cross(A,B);
    auto numerator1 = m + cross(E, P);
    double d = norm(numerator1) / norm(E);
    if (d <= radius) {
        auto tmp = m + cross(E, P);
        auto numerator = cross(E, tmp);
        auto denominator = norm(E) * norm(E);
        Q = P + 1/denominator * numerator;

        auto wA = norm(cross(Q, B)) / norm(m);
        auto wB = norm(cross(Q, A)) / norm(m);

        bool inside = (wA >= 0) && (wA <= 1) && (wB >= 0) && (wB <= 1) && (d <= radius);
        return inside;
    }

    return false;

}

void check_for_contours(int *neighbours, int *contours, int resolution_x, int resolution_y) {
    int cell_x, cell_y, neighbour;
#pragma omp parallel for private(cell_x, cell_y, neighbour)
    for (cell_x = 1; cell_x < resolution_x - 1; cell_x++) {
        for (cell_y = 1; cell_y < resolution_y - 1; cell_y++) {
            neighbour = neighbours[cell_x * resolution_y + cell_y];

            // Check through all neighbours to see if countour can be set
            if (neighbours[(cell_x - 1) * resolution_y + cell_y - 1] != neighbour ||
                neighbours[cell_x * resolution_y + cell_y - 1] != neighbour ||
                neighbours[(cell_x + 1) * resolution_y + cell_y - 1] != neighbour ||
                neighbours[(cell_x - 1) * resolution_y + cell_y] != neighbour ||
                neighbours[(cell_x + 1) * resolution_y + cell_y] != neighbour ||
                neighbours[(cell_x - 1) * resolution_y + cell_y + 1] != neighbour ||
                neighbours[cell_x * resolution_y + cell_y + 1] != neighbour ||
                neighbours[(cell_x + 1) * resolution_y + cell_y + 1] != neighbour) {
                contours[cell_x * resolution_y + cell_y] = 1;
            } else {
                contours[cell_x * resolution_y + cell_y] = 0;
            }
        }
    }
}


PyArrayObject *make_radial_implementation(PyArrayObject *pos,
                                          PyArrayObject *quant,
                                          PyObject *a,
                                          PyObject *b,
                                          double cylinder_radius,
                                          int nshells) {

    // Guards against segfaults
    if(PyArray_API == nullptr) {
        Py_Initialize();
        import_array()
    }

    npy_intp dims[2] = {2, nshells};
    auto pyProfile = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    int npart   = PyArray_DIMS(pos)[0];
    double a_x  = PyFloat_AsDouble(PyList_GetItem(a, 0));
    double a_y  = PyFloat_AsDouble(PyList_GetItem(a, 1));
    double a_z  = PyFloat_AsDouble(PyList_GetItem(a, 2));
    double b_x  = PyFloat_AsDouble(PyList_GetItem(b, 0));
    double b_y  = PyFloat_AsDouble(PyList_GetItem(b, 1));
    double b_z  = PyFloat_AsDouble(PyList_GetItem(b, 2));

    // Start the clock
    clock_t start = clock();

    auto count = (int *) malloc(nshells * sizeof(int));
    auto profile = (double *) PyArray_DATA(pyProfile);
    memset(profile, 0, 2 * nshells * sizeof(double));
    memset(count, 0, nshells * sizeof(int));

    // A and B denote the 3-vector endpoints of the cylinder axis
    arma::vec A = {a_x, a_y, a_z};
    arma::vec B = {b_x, b_y, b_z};
    double d_total = arma::norm(B - A);

    arma::vec Q;
    arma::vec P;
    double px, py, pz, d_i, d;
    bool inside;

    auto s_pos      = PyArray_STRIDES(pos);
    auto s_quant    = PyArray_STRIDES(quant);

    for (int part = 0; part < npart; part++) {

        px = *(double *) ((char *) PyArray_DATA(pos) + part * s_pos[0] + 0 * s_pos[1]);
        py = *(double *) ((char *) PyArray_DATA(pos) + part * s_pos[0] + 1 * s_pos[1]);
        pz = *(double *) ((char *) PyArray_DATA(pos) + part * s_pos[0] + 2 * s_pos[1]);
        d  = *(double *) ((char *) PyArray_DATA(quant) + part * s_quant[0]);

        P = {px, py, pz};

        inside = inside_cylinder(A, B, P, cylinder_radius, Q);

        // Add a contribution for every point within the cylinder
        if (inside) {
            d_i = arma::norm(Q - A);

            int shell = floor(d_i / d_total * nshells);
            profile[shell]  += d;
            count[shell]    += 1;
        }
    }

    for (int shell = 0; shell < nshells; shell++) {

        // Record the radius of the shell
        profile[nshells + shell] = (double) shell / (double) nshells * d_total;

        // Divide the sum by the count to obtain the average
        if (count[shell] > 0) {
            profile[shell] = profile[shell] / count[shell];
        }
    }

    free(count);
    double time_elapsed = ((double) clock() - (double) start) / CLOCKS_PER_SEC;
    printf("Calculation took %gs\n", time_elapsed);

    return pyProfile;

}

PyObject *make_pcolor_implementation(PyArrayObject *pos,
                                     PyArrayObject *quant,
                                     PyObject *axes,
                                     PyObject *boxsizes,
                                     PyObject *resolutions,
                                     PyObject *centers,
                                     int include_neighbours_in_output,
                                     int numthreads) {

    // Guards against segfaults
    if (PyArray_API == nullptr) {
        Py_Initialize();
        import_array()
    }

    // Pre-set the neighbour and contour containers
    PyArrayObject *pyNeighbours = nullptr;
    PyArrayObject *pyContours = nullptr;
    int *neighbours, *contours;

    int resolution_x = (int) PyLong_AsLong(PyList_GetItem(resolutions, 0));
    int resolution_y = (int) PyLong_AsLong(PyList_GetItem(resolutions, 1));
    double boxsize_x = PyFloat_AsDouble(PyList_GetItem(boxsizes, 0));
    double boxsize_y = PyFloat_AsDouble(PyList_GetItem(boxsizes, 1));
    double boxsize_z = 0.0;
    int axis0 = (int) PyLong_AsLong(PyList_GetItem(axes, 0));
    int axis1 = (int) PyLong_AsLong(PyList_GetItem(axes, 1));
    int axis2 = 3 - axis0 - axis1;
    double center_x = PyFloat_AsDouble(PyList_GetItem(centers, 0));
    double center_y = PyFloat_AsDouble(PyList_GetItem(centers, 1));
    double center_z = PyFloat_AsDouble(PyList_GetItem(centers, 2));

    // Initialise other objects
    npy_intp dims[2] = {resolution_x, resolution_y};

    int npart = PyArray_DIMS(pos)[0];
    auto pyGrid = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    auto grid = (double *) PyArray_DATA(pyGrid);

    // Initialise grid to zero
    memset(grid, 0, resolution_x * resolution_y * sizeof(double));
    double cellsizex = boxsize_x / resolution_x;
    double cellsizey = boxsize_y / resolution_y;
    double cellsizez = 0;

    // Start the timer
    double start = get_time();


    // Initialise array objects for neighbours and countours
    if (include_neighbours_in_output) {
        pyNeighbours = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_INT);
        pyContours = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_INT);
        neighbours = (int *) PyArray_DATA(pyNeighbours);
        contours = (int *) PyArray_DATA(pyContours);
        memset(contours, 0, resolution_x * resolution_y * sizeof(int));
    }

    // Create a local container for position and quantity
    auto real_pos = (double *) malloc(3 * npart * sizeof(double));
    auto real_quant = (double *) malloc(npart * sizeof(double));


    int i, j;
#pragma omp parallel for private(i, j)
    for (i = 0; i < npart; i++) {
        for (j = 0; j < 3; j++) {
            real_pos[i * 3 + j] = *(double *) ((char *) PyArray_DATA(pos) + i * PyArray_STRIDES(pos)[0] +
                                               j * PyArray_STRIDES(pos)[1]);
        }
        real_quant[i] = *(double *) ((char *) PyArray_DATA(quant) + i * PyArray_STRIDES(quant)[0]);
    }

    // Create a tree using the position info
    t_sph_tree tree;
    createTree(&tree, npart, real_pos);

    printf("Starting tree walk with %d thread(s)\n", numthreads);
    set_num_threads(numthreads);

    double start_tree = get_time();

    int x, y, z;
    int neighbour = -1;
    int cell = 0;
#pragma omp parallel private(x, y, z) firstprivate(cell, neighbour)
    {
        int *worklist = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));

        double coord[3];
        int count = 0;
        int thread_id = get_thread_id();

        coord[axis2] = center_z;

#pragma omp for schedule(dynamic, 1) nowait
        for (x = 0; x < resolution_x; x++) {
            coord[axis0] = center_x - 0.5 * boxsize_x + cellsizex * (0.5 + x);
            for (y = 0; y < resolution_y; y++) {
                cell = y + resolution_y * x;
                coord[axis1] = center_y - 0.5 * boxsize_y + cellsizey * (0.5 + y);
                for (z = 0; z < 1; z++) {
                    coord[axis2] = center_z - 0.5 * boxsize_z + cellsizez * (0.5 + z);

                    getNearestNeighbour(&tree, real_pos, coord, &neighbour, worklist);

                    //
#pragma omp atomic
                    grid[cell] += real_quant[neighbour];

                    if (include_neighbours_in_output) {
                        neighbours[cell] = neighbour;
                    }
                }
            }

            if (include_neighbours_in_output) {
                neighbour = neighbours[cell - resolution_y];
            }

            if (thread_id == numthreads - 1) {
                if (x / (resolution_x / 10) > count) {
                    count++;
                    printf("Done iter %4d of %4d\n", x + 1, resolution_x);
                }
            }
        }

        free(worklist);
    }

    // Write info about tree walk
    double tree_walk_time = get_time() - start_tree;
    printf("Tree walk took %gs\n", tree_walk_time);

    // Free memory
    freeTree(&tree);
    free(real_pos);
    free(real_quant);

    // Create contours using the neighbour information
    if (include_neighbours_in_output) {
        check_for_contours(neighbours, contours, resolution_x, resolution_y);
    }

    // Write to dict
    PyObject *dict = PyDict_New();
    PyDict_SetStolenItem(dict, "grid", (PyObject *) pyGrid);
    if (include_neighbours_in_output) {
        PyDict_SetStolenItem(dict, "neighbours", (PyObject *) pyNeighbours);
        PyDict_SetStolenItem(dict, "contours", (PyObject *) pyContours);
    }

    // Show the elapsed time
    double calculation_time = get_time() - start;
    printf("Calculation took %gs\n", calculation_time);

    return dict;
}

PyObject *make_pcolor(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *pos, *quant;
    PyObject *axes, *boxsizes, *resolutions, *centers;
    int numthreads = 1;
    int include_neighbours_in_output = 1; // Don't try to turn the include variable into a bool. It caused nothing but problems

    const char *kwlist[] = {"pos", "quant", "axes", "boxsizes", "resolutions", "centers", "include_neighbours_in_output", "numthreads", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O!O!OOOO|ii:make_pcolor("
                                     "pos, quant, "
                                     "axes, "
                                     "boxsizes, "
                                     "resolutions, "
                                     "centers, "
                                     "[include_neighbours_in_output, "
                                     "numthreads])",
                                     keywords,
                                     &PyArray_Type, &pos,
                                     &PyArray_Type, &quant,
                                     &axes,
                                     &boxsizes,
                                     &resolutions,
                                     &centers,
                                     &include_neighbours_in_output,
                                     &numthreads)) {
        return nullptr;
    }


    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    if (PyArray_NDIM(quant) != 1 || PyArray_TYPE(quant) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "quant has to be of dimension [n] and type double");
        return nullptr;
    }


    PyObject *dict = make_pcolor_implementation(pos, quant, axes, boxsizes, resolutions, centers,
                                                include_neighbours_in_output, numthreads);

    return dict;
}

PyObject *make_radial(PyObject *self, PyObject *args) {

    PyArrayObject *pos, *quant;
    PyObject *a, *b;
    double cylinder_radius;
    int nshells = 0;

    if (!PyArg_ParseTuple(args, "O!O!OOdi:make_radial(pos, quant, a, b, cylinder_radius, nshells)",
                          &PyArray_Type, &pos, &PyArray_Type, &quant,
                          &a, &b,
                          &cylinder_radius,
                          &nshells)) {
        return nullptr;
    }

    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    if (PyArray_NDIM(quant) != 1 || PyArray_TYPE(quant) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "quant has to be of dimension [n] and type double");
        return nullptr;
    }

    auto pyProfile = make_radial_implementation(pos, quant, a, b, cylinder_radius, nshells);

    return PyArray_Return(pyProfile);
}

PyObject *get_indices_of_neighbours(PyObject *self, PyObject *args) {

    PyObject *dict;
    int n_neighbours;
    PyArrayObject *target_coord;

    if (!PyArg_ParseTuple(args, "OiO!:get_indices_of_neighbours(dict, n_neighbours, target_coord)",
                          &dict,
                          &n_neighbours,
                          &PyArray_Type, &target_coord)) {
        return nullptr;
    }

    auto pos = (double *) PyArray_DATA((PyArrayObject *) PyDict_GetItemString(dict, f[N::COORDINATES]));
    auto coord = (double *) PyArray_DATA(target_coord);
    auto mass = (double *) PyArray_DATA((PyArrayObject *) PyDict_GetItemString(dict, f[N::MASSES]));
    auto density_1 = (PyArrayObject *) PyDict_GetItemString(dict, f[N::DENSITY]);
    auto density = (double *) PyArray_DATA(density_1);
    int npart = PyArray_DIMS(density_1)[0];

    t_sph_tree tree;
    createTree(&tree, npart, pos);

    double hsml = 0;
    calcHsml(&tree, coord, pos, mass, n_neighbours, &hsml, density);

    int nneighbours_real;
    int converged;
    auto neighbours = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));

    auto radius = getNNeighbours(&tree, coord, pos,
                                 n_neighbours,
                                 &nneighbours_real,
                                 &neighbours,
                                 &converged);

    if (nneighbours_real == 0 && converged == 1) {
        coord[0] += 0.0001 * coord[0];
        coord[1] += 0.0001 * coord[1];
        coord[2] += 0.0001 * coord[2];

        radius = getNNeighbours(&tree, coord, pos,
                                n_neighbours,
                                &nneighbours_real,
                                &neighbours,
                                &converged);
    }

    auto ret_dict = PyDict_New();
    PyDict_SetStolenItem(ret_dict, "Radius", (PyObject*) PyFloat_FromDouble(radius));
    PyDict_SetStolenItem(ret_dict, "Converged", (PyObject *) PyLong_FromLong(converged));
    PyDict_SetStolenItem(ret_dict, "NNeighboursReal", (PyObject *) PyLong_FromLong(nneighbours_real));
    PyDict_SetStolenItem(ret_dict, "Neighbours", (PyObject *) create_numpy_array(neighbours, nneighbours_real));

    freeTree(&tree);
    free(neighbours);

    return ret_dict;
}


PyMODINIT_FUNC PyInit_arepo_vis(void) {
    import_array()

    return PyModule_Create(&moduledef_arepo_vis);
}