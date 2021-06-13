#include <Python.h>
#include <arrayobject.h>

#include "sph.h"
#include "utils.h"
#include "omp_util.h"
#include "make_pcolor.h"



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

    PyArrayObject *pos, *quant, *axes, *boxsizes, *resolutions, *centers;
    int numthreads = 1;
    int include_neighbours_in_output = 1; // Don't try to turn the include variable into a bool. It caused nothing but problems

    const char *kwlist[] = {"pos", "quant", "axes", "boxsizes", "resolutions", "centers",
                            "include_neighbours_in_output", "numthreads", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O!O!O!O!O!O!|ii:make_pcolor("
                                     "pos, quant, "
                                     "axes, "
                                     "boxsizes, "
                                     "resolutions, "
                                     "centers, "
                                     "[include_neighbours_in_output, "
                                     "numthreads] )",
                                     keywords,
                                     &PyArray_Type, &pos, &PyArray_Type, &quant,
                                     &PyArray_Type, &axes,
                                     &PyArray_Type, &boxsizes,
                                     &PyArray_Type, &resolutions,
                                     &PyArray_Type, &centers,
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

PyObject *
make_pcolor_implementation(PyArrayObject *pos,
                           PyArrayObject *quant,
                           PyArrayObject *axes,
                           PyArrayObject *boxsizes,
                           PyArrayObject *resolutions,
                           PyArrayObject *centers,
                           int include_neighbours_in_output,
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

    int resolution_x    = * (int *) ((char *) PyArray_DATA(resolutions) + 0 * PyArray_STRIDES(resolutions)[0]);
    int resolution_y    = * (int *) ((char *) PyArray_DATA(resolutions) + 1 * PyArray_STRIDES(resolutions)[0]);
    double boxsize_x    = * (double *) ((char *) PyArray_DATA(boxsizes) + 0 * PyArray_STRIDES(boxsizes)[0]);
    double boxsize_y    = * (double *) ((char *) PyArray_DATA(boxsizes) + 1 * PyArray_STRIDES(boxsizes)[0]);
    double boxsize_z    = 0;
    int axis0           = * (int *) ((char *) PyArray_DATA(axes) + 0 * PyArray_STRIDES(axes)[0]);
    int axis1           = * (int *) ((char *) PyArray_DATA(axes) + 1 * PyArray_STRIDES(axes)[0]);
    int axis2           = 3 - axis0 - axis1;
    double center_x     = * (double *) ((char *) PyArray_DATA(centers) + axis0 * PyArray_STRIDES(centers)[0]);
    double center_y     = * (double *) ((char *) PyArray_DATA(centers) + axis1 * PyArray_STRIDES(centers)[0]);
    double center_z     = * (double *) ((char *) PyArray_DATA(centers) + axis2 * PyArray_STRIDES(centers)[0]);

    // Initialise other objects
    npy_intp dims[2] = {resolution_x, resolution_y};

    int npart           = PyArray_DIMS(pos)[0];
    auto pyGrid         = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    auto grid           = (double *) PyArray_DATA(pyGrid);

    // Initialise grid to zero
    memset(grid, 0, resolution_x * resolution_y * sizeof(double));
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
        real_quant[i] = *(double *) ((char *) PyArray_DATA(quant) + i * PyArray_STRIDES(quant)[0]);
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

// Python Module definition
static PyMethodDef pcolor_methods[] = {
        {"make_pcolor", (PyCFunction) make_pcolor, METH_VARARGS | METH_KEYWORDS, ""},
        {nullptr, nullptr, 0,nullptr}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "arepo_pcolor",
        nullptr,
        -1,
        pcolor_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_arepo_pcolor(void) {
    import_array();

    return PyModule_Create(&moduledef);
}
