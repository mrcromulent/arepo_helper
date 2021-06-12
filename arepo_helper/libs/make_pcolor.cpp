#include <Python.h>
#include <arrayobject.h>

#include "sph.h"
#include "omp_util.h"
#include "make_pcolor.h"
#include "utils.h"


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


PyObject *make_pcolor(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *pos, *value;
    int resolution_x, resolution_y;
    double boxsize_x, boxsize_y, boxsize_z;
    double center_x, center_y, center_z;
    int axis0 = 0, axis1 = 1;
    int nz = 1, numthreads = 1;
    int include_neighbours_in_output = 1;

    char *kwlist[] = {"pos", "value",
                      "resolution_x", "resolution_y",
                      "boxsize_x", "boxsize_y", "boxsize_z",
                      "center_x", "center_y", "center_z",
                      "axis0", "axis1",
                      "nz",
                      "include_neighbours_in_output",
                      "numthreads", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O!O!iidddddd|iiiii:make_pcolor("
                                     "pos, value, "
                                     "resolution_x, resolution_y, "
                                     "boxsize_x, boxsize_y, boxsize_z, "
                                     "centerx, centery, centerz, "
                                     "[axis0, axis1, "
                                     "nz, "
                                     "include_neighbours_in_output, "
                                     "numthreads"
                                     "] )",
                                     kwlist,
                                     &PyArray_Type, &pos, &PyArray_Type, &value,
                                     &resolution_x, &resolution_y,
                                     &boxsize_x, &boxsize_y, &boxsize_z,
                                     &center_x, &center_y, &center_z,
                                     &axis0, &axis1,
                                     &nz,
                                     &include_neighbours_in_output,
                                     &numthreads)) {
        return nullptr;
    }


    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "value has to be of dimension [n] and type double");
        return nullptr;
    }

//    std::cout << "resolution_x = " << resolution_x << std::endl;


    PyObject *dict = make_pcolor_implementation(pos, value,
                                           resolution_x, resolution_y,
                                           boxsize_x, boxsize_y, boxsize_z,
                                           center_x, center_y, center_z,
                                           axis0, axis1,
                                           nz,
                                           include_neighbours_in_output,
                                           numthreads);

//    auto grid = (PyArrayObject *) PyDict_GetItemString(dict, "grid");
//    auto grid_data = (double *) PyArray_DATA(grid);
//    auto n = PyArray_DIMS(grid)[0];

//    std::cout << n << std::endl;
//    for (int i = 0; i < n; i++) {
//        if (i % 10 == 0) {
//            std::cout << grid_data[i] << std::endl;
//        }
//    }

    return dict;
}

PyObject *make_pcolor_implementation(PyArrayObject *pos, PyArrayObject *value,
                                     int resolution_x, int resolution_y,
                                     double boxsize_x, double boxsize_y, double boxsize_z,
                                     double center_x, double center_y, double center_z,
                                     int axis0, int axis1,
                                     int nz,
                                     bool include_neighbours_in_output,
                                     int numthreads) {

    // Guards against segfaults
    if(PyArray_API == nullptr) {
        Py_Initialize();
        import_array()
    }

    // Pre-set the neighbour and contour containers
    PyArrayObject *pyNeighbours = nullptr;
    PyArrayObject *pyContours = nullptr;
    int *neighbours, *contours;

    // Initialise other objects
    npy_intp dims[2] = {resolution_x, resolution_y};
    int axis2           = 3 - axis0 - axis1;
    int npart           = PyArray_DIMS(pos)[0];
    auto pyGrid         = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    auto grid           = (double *) PyArray_DATA(pyGrid);

    // Initialise grid to zero
    memset(grid, 0, resolution_x * resolution_y * sizeof(double));
//    std::cout << "resolution_x = " << resolution_x << std::endl;

    double cellsizex    = boxsize_x / resolution_x;
    double cellsizey    = boxsize_y / resolution_y;
    double cellsizez    = 0;

    // Start the timer
    double start = get_time();


    // Initialise array objects for neighbours and countours
    if (include_neighbours_in_output) {
        pyNeighbours    = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_INT);
        pyContours      = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_INT);
        neighbours      = (int *) PyArray_DATA(pyNeighbours);
        contours        = (int *) PyArray_DATA(pyContours);
        memset(contours, 0, resolution_x * resolution_y * sizeof(int));
    }

    // Create a local container for position and quantity
    auto real_pos = (double *) malloc(3 * npart * sizeof(double));
    auto real_quant = (double *) malloc(npart * sizeof(double));


    int i,j;
#pragma omp parallel for private(i, j)
    for (i = 0; i < npart; i++) {
        for (j = 0; j < 3; j++) {
            real_pos[i * 3 + j] = *(double *) ((char *) PyArray_DATA(pos) + i * PyArray_STRIDES(pos)[0] + j * PyArray_STRIDES(pos)[1]);
        }
        real_quant[i] = *(double *) ((char *) PyArray_DATA(value) + i * PyArray_STRIDES(value)[0]);
    }

    // Create a tree using the position info
    t_sph_tree tree;
    createTree(&tree, npart, real_pos);

    printf("Starting tree walk with %d thread(s)\n", numthreads);
    set_num_threads(numthreads);

    double start_tree = get_time();

    int x, y, z;
    int neighbour   = -1;
    int cell        = 0;
#pragma omp parallel private(x, y, z) firstprivate(cell, neighbour)
    {
        int *worklist = static_cast<int *>(malloc(tree.usednodes * sizeof(int)));

        double coord[3];
        int count = 0;
        int thread_id = get_thread_id();

        coord[axis2] = center_z;

#pragma omp for schedule(dynamic, 1) nowait
        for (x = 0; x < resolution_x; x++) {
//            std::cout << "x = " << x << std::endl;
            coord[axis0] = center_x - 0.5 * boxsize_x + cellsizex * (0.5 + x);
            for (y = 0; y < resolution_y; y++) {
                cell = y + resolution_y * x;
                coord[axis1] = center_y - 0.5 * boxsize_y + cellsizey * (0.5 + y);
                for (z = 0; z < nz; z++) {
                    coord[axis2] = center_z - 0.5 * boxsize_z + cellsizez * (0.5 + z);

                    getNearestNeighbour(&tree, real_pos, coord, &neighbour, worklist);

                    //
#pragma omp atomic
                    grid[cell] += real_quant[neighbour];
//                    std::cout << real_quant[neighbour] << std::endl;

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

    auto max_val = resolution_x * resolution_y;
//    std::cout << "Just after tree walk" << std::endl;
//    std::cout << grid[0] << std::endl;
//    for (int i = 0; i < max_val; i++) {
//        if (i % 10000 == 0) {
//            std::cout << "i = " << i << ". q = " << grid[i] << std::endl;
//        }
//    }

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
    auto success = PyDict_SetStolenItem(dict, "grid", (PyObject *) pyGrid);
//    PyDict_SetItemString(dict, "grid2", (PyObject *) 1);
    if (include_neighbours_in_output) {
        PyDict_SetStolenItem(dict, "neighbours", (PyObject *) pyNeighbours);
        PyDict_SetStolenItem(dict, "contours", (PyObject *) pyContours);
    }

    // Show the elapsed time
    double calculation_time = get_time() - start;
    printf("Calculation took %gs\n", calculation_time);

    auto grid2 = (PyArrayObject *) PyDict_GetItemString(dict, "grid");
    auto grid_data = (double *) PyArray_DATA(grid2);
    auto n = PyArray_DIMS(grid2)[0];

//    std::cout << n << std::endl;
//    for (int i = 0; i < n; i++) {
//        if (i % 10 == 0) {
//            std::cout << grid_data[i] << std::endl;
//        }
//    }

    return dict;
}

// Python Module definition
static PyMethodDef pcolor_methods[] = {
        {"make_pcolor", (PyCFunction) make_pcolor, METH_VARARGS | METH_KEYWORDS, ""},
        {nullptr, nullptr, 0,nullptr}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pcolor_pierre",
        nullptr,
        -1,
        pcolor_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_pcolor_pierre(void) {
    import_array();

    return PyModule_Create(&moduledef);
}
