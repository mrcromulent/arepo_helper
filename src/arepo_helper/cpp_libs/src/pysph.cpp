#include <Python.h>
#include <arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "omp_util.h"
#include "pysph.h"

inline void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
    PyDict_SetItemString(dict, key, object);
    Py_DECREF(object);
}

PyArrayObject* createPyArray2Dim(double *data, int dim1, int dim2) {
    PyArrayObject* pyData;

    if (dim1 == 1 || dim2 == 1) {
        npy_intp dims[1];
        dims[0] = dim1 * dim2;
        pyData = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        memcpy(PyArray_DATA(pyData), data, dim1*dim2*sizeof(double));
    } else {
        npy_intp dims[2];
        dims[0] = dim1;
        dims[1] = dim2;
        pyData = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        memcpy(PyArray_DATA(pyData), data, dim1*dim2*sizeof(double));
    }

    return pyData;
}

PyArrayObject* createPyArray(double *data, int length) {
    PyArrayObject* pyData;

    npy_intp dims[1];
    dims[0] = length;

    pyData = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    memcpy(PyArray_DATA(pyData), data, length*sizeof(double));

    return pyData;
}

PyArrayObject* createPyIntArray(int *data, int length) {
    PyArrayObject* pyData;

    npy_intp dims[1];
    dims[0] = length;

    pyData = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    memcpy(PyArray_DATA(pyData), data, length*sizeof(int));

    return pyData;
}

typedef struct {
    PyObject_HEAD;
    t_sph_tree *tree;
} pyTree;

/* forward declarations */
void pyTreeDealloc(PyObject* self);
static PyObject* pyTree_str(pyTree* self);
PyObject* pyCalcDensity(pyTree* self, PyObject* args);
PyObject* pyCalcHsml(pyTree* self, PyObject* args);
PyObject* pyCalcHsmlMulti(pyTree* self, PyObject* args, PyObject* kwargs);
PyObject* pyGetNearestNeighbours(pyTree* self, PyObject* args);

static PyMethodDef pyTreeMethods[] =
        {
                { "calcDensity", (getattrofunc)pyCalcDensity, METH_VARARGS, nullptr },
                { "calcHsml", (getattrofunc)pyCalcHsml, METH_VARARGS, nullptr },
                { "calcHsmlMulti", (getattrofunc)pyCalcHsmlMulti, METH_VARARGS | METH_KEYWORDS, nullptr },
                { "getNearestNeighbours", (getattrofunc)pyGetNearestNeighbours, METH_VARARGS, nullptr },
                { nullptr, nullptr, 0, nullptr }
        };


static PyTypeObject pyTreeType = {
        PyVarObject_HEAD_INIT(nullptr, 0)
        "pyeos",               /* tp_name */
        sizeof(pyTree),         /* tp_basicsize */
        0,                         /* tp_itemsize */
        pyTreeDealloc,                         /* tp_dealloc */
        0,                         /* tp_print */
        nullptr,                         /* tp_getattr */
        nullptr,                         /* tp_setattr */
        nullptr,                         /* tp_reserved */
        nullptr,                         /* tp_repr */
        nullptr,                         /* tp_as_number */
        nullptr,                         /* tp_as_sequence */
        nullptr,                         /* tp_as_mapping */
        nullptr,                         /* tp_hash  */
        nullptr,                         /* tp_call */
        (reprfunc)pyTree_str,                         /* tp_str */
        nullptr,                         /* tp_getattro */
        nullptr,                         /* tp_setattro */
        nullptr,                         /* tp_as_buffer */
        Py_TPFLAGS_DEFAULT,        /* tp_flags */
        "pytree",            /* tp_doc */
        nullptr,                         /* tp_traverse */
        nullptr,                         /* tp_clear */
        nullptr,                         /* tp_richcompare */
        0,                         /* tp_weaklistoffset */
        nullptr,                         /* tp_iter */
        nullptr,                         /* tp_iternext */
        pyTreeMethods,          /* tp_methods */
        nullptr,                         /* tp_members */
        nullptr,                         /* tp_getset */
        nullptr,                         /* tp_base */
        nullptr,                         /* tp_dict */
        nullptr,                         /* tp_descr_get */
        nullptr,                         /* tp_descr_set */
        0,                         /* tp_dictoffset */
        nullptr,                         /* tp_init */
        nullptr,                         /* tp_alloc */
        nullptr,                         /* tp_new */
};

void pyTreeDealloc(PyObject* self) {
    freeTree(((pyTree*)self)->tree);
    PyObject_Del(self);
}

PyObject* pyTree_str(pyTree* self) {
    return PyUnicode_FromString("tree");
}

PyObject* pyCalcDensity(pyTree* self, PyObject* args) {
    PyArrayObject *pyCoord, *pyPos, *pyMass;
    double *coord, *pos, *mass, *real_pos;
    double hsml, density;
    int npart, i, j;

    if (!PyArg_ParseTuple(args, "O!dO!O!:calcDensity(coord, hsml, pos, mass)", &PyArray_Type, &pyCoord, &hsml, &PyArray_Type, &pyPos, &PyArray_Type, &pyMass)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyCoord) != 1 || PyArray_DIMS(pyCoord)[0] != 3 || PyArray_TYPE(pyCoord) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "coord has to be of dimensions [3] and type double");
        return nullptr;
    }
    coord = (double*)PyArray_DATA(pyCoord);

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }
    pos = (double*)PyArray_DATA(pyPos);
    npart = PyArray_DIMS(pyPos)[0];

    if (PyArray_NDIM(pyMass) != 1 || PyArray_DIMS(pyMass)[0] != npart || PyArray_TYPE(pyMass) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "mass has to be of dimensions [n] and type double");
        return nullptr;
    }
    mass = (double*)PyArray_DATA(pyMass);

    if (npart != self->tree->npart) {
        PyErr_SetString(PyExc_ValueError, "tree was calculated with the wrong number of particles.");
        return nullptr;
    }

    real_pos = (double*)malloc(3 * npart * sizeof(double));
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)pos + i*PyArray_STRIDES(pyPos)[0] + j*PyArray_STRIDES(pyPos)[1]);
        }
    }

    calcDensity(self->tree, coord, hsml, real_pos, mass, &density, nullptr);

    free(real_pos);

    return Py_BuildValue("d", density);
}

PyObject* pyCalcHsml(pyTree* self, PyObject* args) {
    PyArrayObject *pyCoord, *pyPos, *pyMass;
    double *mass, *real_pos;
    double center[3];
    double hsml, density, weighted_neighbours;
    int npart, nneighbours, i, j;

    hsml = 0;
    if (!PyArg_ParseTuple(args, "O!O!O!i|d:calcHsml(coord, pos, mass, neighbours, [hsml])", &PyArray_Type, &pyCoord, &PyArray_Type, &pyPos, &PyArray_Type, &pyMass, &nneighbours, &hsml)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyCoord) != 1 || PyArray_DIMS(pyCoord)[0] != 3 || PyArray_TYPE(pyCoord) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "coord has to be of dimensions [3] and type double");
        return nullptr;
    }
    for (i=0; i<3; i++) {
        center[i] = *(double*)((char*)PyArray_DATA(pyCoord) + i*PyArray_STRIDES(pyCoord)[0]);
    }

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    npart = PyArray_DIMS(pyPos)[0];
    if (PyArray_NDIM(pyMass) != 1 || PyArray_DIMS(pyMass)[0] != npart || PyArray_TYPE(pyMass) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "mass has to be of dimensions [n] and type double");
        return nullptr;
    }
    mass = (double*)PyArray_DATA(pyMass);

    if (npart != self->tree->npart) {
        PyErr_SetString(PyExc_ValueError, "tree was calculated with the wrong number of particles.");
        return nullptr;
    }

    real_pos = (double*)malloc(3 * npart * sizeof(double));
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pyPos) + i*PyArray_STRIDES(pyPos)[0] + j*PyArray_STRIDES(pyPos)[1]);
        }
    }

    weighted_neighbours = calcHsml(self->tree, center, real_pos, mass, nneighbours, &hsml, &density);

    free(real_pos);

    return Py_BuildValue("ddd", hsml, density, weighted_neighbours);
}

PyObject* pyCalcHsmlMulti(pyTree* self, PyObject* args, PyObject* kwargs) {
    PyArrayObject *pyCoord, *pyPos, *pyMass, *pyHsml, *pyHsmlGuess;
    double *mass, *hsml, *real_pos, *hsmlGuess;
    double density;
    int npart, nneighbours, multi, i, j, numthreads;
    double start;
    const char *kwlist[] = {"coord", "pos", "mass", "neighbours", "hsmlGuess", "numthreads", nullptr};
    auto keywords = (char **) kwlist;

    start = get_time();

    pyHsmlGuess = nullptr;
    numthreads = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!i|O!i:calcHsmlMulti(coord, pos, mass, neighbours, [hsmlGuess, numthreads])", keywords, &PyArray_Type, &pyCoord, &PyArray_Type, &pyPos, &PyArray_Type, &pyMass, &nneighbours, &PyArray_Type, &pyHsmlGuess, &numthreads)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyCoord) != 2 || PyArray_DIMS(pyCoord)[1] != 3 || PyArray_TYPE(pyCoord) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "coord has to be of dimensions [m,3] and type double");
        return nullptr;
    }
    multi = PyArray_DIMS(pyCoord)[0];

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }
    npart = PyArray_DIMS(pyPos)[0];

    if (PyArray_NDIM(pyMass) != 1 || PyArray_DIMS(pyMass)[0] != npart || PyArray_TYPE(pyMass) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "mass has to be of dimensions [n] and type double");
        return nullptr;
    }
    mass = (double*)PyArray_DATA(pyMass);

    hsmlGuess = nullptr;
    if (pyHsmlGuess) {
        hsmlGuess = (double*)PyArray_DATA(pyHsmlGuess);
        if (PyArray_NDIM(pyHsmlGuess) != 1 || PyArray_DIMS(pyHsmlGuess)[0] != multi || PyArray_TYPE(pyHsmlGuess) != NPY_DOUBLE) {
            PyErr_SetString(PyExc_ValueError, "hsmlGuess has to be of dimensions [m] and type double");
            return nullptr;
        }
    }

    if (npart != self->tree->npart) {
        PyErr_SetString(PyExc_ValueError, "tree was calculated with the wrong number of particles.");
        return nullptr;
    }

    set_num_threads(numthreads);

    real_pos = (double*)malloc(3 * npart * sizeof(double));

#pragma omp parallel for private(i, j)
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pyPos) + i*PyArray_STRIDES(pyPos)[0] + j*PyArray_STRIDES(pyPos)[1]);
        }
    }

    printf("Doing %d positions with %d thread(s).\n", multi, numthreads);

    hsml = (double*)malloc(multi * sizeof(double));
    if (pyHsmlGuess)
    {
        printf("HsmlGuess found.\n");
        memcpy(hsml, hsmlGuess, multi * sizeof(double));
    }
    else
    {
        printf("No hsmlGuess found.\n");
        memset(hsml, 0, multi * sizeof(double));
    }

#pragma omp parallel private(i, j, density)
    {
        int slowly = 0;
        int nextoutput = multi / 100;
        double runtime, center[3];
        int thread_id = get_thread_id();

#pragma omp for schedule(dynamic, 1) nowait
        for (i=0; i<multi; i++) {
            for (j=0; j<3; j++) {
                center[j] = *(double*)((char*)PyArray_DATA(pyCoord) + i*PyArray_STRIDES(pyCoord)[0] + j*PyArray_STRIDES(pyCoord)[1]);
            }

            calcHsml(self->tree, center, real_pos, mass, nneighbours, &hsml[i], &density);

            if (thread_id == numthreads - 1){
                if (i >= nextoutput) {
                    runtime = get_time() - start;

                    if (nextoutput == multi / 100) {
                        if (runtime > 60.) {
                            slowly = 1;
                        } else {
                            nextoutput = 0;
                        }
                    }

                    printf("%d / %d particles done (%d%%): %ds elapsed, ~%ds remaining\n", i, multi, (int)floor(100.0*(double)i/(double)multi), (int)(runtime), (int)(runtime/i*(multi-i)));

                    if (slowly)
                        nextoutput += multi / 100;
                    else
                        nextoutput += multi /  10;
                }
            }
        }
    }

    free(real_pos);

    printf("Hsml calculation took %gs\n", get_time() - start);

    pyHsml = createPyArray(hsml, multi);
    free(hsml);

    return (PyObject*)pyHsml;
}


PyObject* pyGetNearestNeighbours(pyTree* self, PyObject* args) {
    PyArrayObject *pyCoord, *pyPos, *pyNeighbours;
    double *pos, *real_pos;
    int npart, multi, i, j, *neighbours;
    int numthreads;
    time_t start;

    start = clock();

    numthreads = 1;
    if (!PyArg_ParseTuple(args, "O!O!|i:pyGetNearestNeighbours(coord, pos, [numthreads])", &PyArray_Type, &pyCoord, &PyArray_Type, &pyPos, &numthreads)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyCoord) != 2 || PyArray_DIMS(pyCoord)[1] != 3 || PyArray_TYPE(pyCoord) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "coord has to be of dimensions [m,3] and type double");
        return nullptr;
    }
    multi = PyArray_DIMS(pyCoord)[0];

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }
    pos = (double*)PyArray_DATA(pyPos);
    npart = PyArray_DIMS(pyPos)[0];

    if (npart != self->tree->npart) {
        PyErr_SetString(PyExc_ValueError, "tree was calculated with the wrong number of particles.");
        return nullptr;
    }

    set_num_threads(numthreads);

    real_pos = (double*)malloc(3 * npart * sizeof(double));

#pragma omp parallel for private(i, j)
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)pos + i*PyArray_STRIDES(pyPos)[0] + j*PyArray_STRIDES(pyPos)[1]);
        }
    }

    neighbours = (int*)malloc(multi * sizeof(int));

#pragma omp parallel private(i, j)
    {
        int *worklist = static_cast<int *>(malloc(self->tree->usednodes * sizeof(int)));

        int slowly = 0;
        int nextoutput = multi / 100;
        double runtime;
        int thread_id = get_thread_id();
        double center[3];

#pragma omp for schedule(dynamic, 1000) nowait
        for (i=0; i<multi; i++) {
            for (j=0; j<3; j++) {
                center[j] = *(double*)((char*)PyArray_DATA(pyCoord) + i*PyArray_STRIDES(pyCoord)[0] + j*PyArray_STRIDES(pyCoord)[1]);
            }
            neighbours[i] = getNearestNeighbour(self->tree, pos, center, nullptr, worklist);

            if (thread_id == numthreads - 1) {
                if (i >= nextoutput) {
                    runtime = ((double)clock()-(double)start)/CLOCKS_PER_SEC;

                    if (nextoutput == multi / 100) {
                        if (runtime > 60.) {
                            slowly = 1;
                        } else {
                            nextoutput = 0;
                        }
                    }

                    printf("%d / %d coordinates done (%d%%): %ds elapsed, ~%ds remaining\n", i, multi, (int)floor(100.0*(double)i/(double)multi), (int)(runtime), (int)(runtime/i*(multi-i)));
                    if (slowly)
                        nextoutput += multi / 100;
                    else
                        nextoutput += multi /  10;
                }
            }
        }

        free(worklist);
    }
    free(real_pos);

    pyNeighbours = createPyIntArray(neighbours, multi);
    free(neighbours);

    printf("Neighbour search took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC);

    return (PyObject*)pyNeighbours;
}

PyObject* PyGetTreeObject(t_sph_tree* tree) {
    pyTree* pytree;

    if (tree) {
        pytree = (pyTree*)PyObject_New(pyTree, &pyTreeType);
        pytree->tree = tree;
        return (PyObject*)pytree;
    } else {
        Py_RETURN_NONE;
    }
}

int pyConvertTree(PyObject* object, t_sph_tree** tree) {
    if (strcmp(object->ob_type->tp_name, "pytree") != 0) {
        PyErr_BadArgument();
        return 0;
    }

    *tree = ((pyTree*)object)->tree;
    return 1;
}

/* pyTree module */
PyObject* makeTree(PyObject *self, PyObject *args) {
    t_sph_tree* tree;
    PyArrayObject *pos;
    double *data_pos, *real_pos;
    int npart, i, j;

    if (!PyArg_ParseTuple(args, "O!:makeTree(pos)", &PyArray_Type, &pos)) {
        return nullptr;
    }

    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }
    data_pos = (double*)PyArray_DATA(pos);
    npart = PyArray_DIMS(pos)[0];

    real_pos = (double*)malloc(3 * npart * sizeof(double));
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)data_pos + i*PyArray_STRIDES(pos)[0] + j*PyArray_STRIDES(pos)[1]);
        }
    }

    tree = (t_sph_tree*)malloc(sizeof(t_sph_tree));
    createTree(tree, npart, real_pos);

    free(real_pos);

    return PyGetTreeObject(tree);
}

PyObject* relaxData(PyObject *self, PyObject *args) {
    PyObject *wd;
    PyArrayObject *pyPos, *pyMass, *pyHsml, *wdRho, *pyNewPos;
    double *mass, *hsml, *real_pos;
    double *rho, dr;
    int nneighbours, nsteps;
    int npart, ncells, i, j;

    if (!PyArg_ParseTuple(args, "O!O!O!iOi:relaxData(pos, mass, hsml, neighbours, wd, nsteps)", &PyArray_Type, &pyPos, &PyArray_Type, &pyMass, &PyArray_Type, &pyHsml, &nneighbours, &wd, &nsteps)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }
    npart = PyArray_DIMS(pyPos)[0];

    real_pos = (double*)malloc(3 * npart * sizeof(double));
    for (i=0; i<npart; i++) {
        for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pyPos) + i*PyArray_STRIDES(pyPos)[0] + j*PyArray_STRIDES(pyPos)[1]);
        }
    }

    if (PyArray_NDIM(pyMass) != 1 || PyArray_DIMS(pyMass)[0] != npart || PyArray_TYPE(pyMass) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "mass has to be of dimensions [n] and type double");
        return nullptr;
    }
    mass = (double*)PyArray_DATA(pyMass);

    if (PyArray_NDIM(pyHsml) != 1 || PyArray_DIMS(pyHsml)[0] != npart || PyArray_TYPE(pyHsml) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "mass has to be of dimensions [n] and type double");
        return nullptr;
    }
    hsml = (double*)PyArray_DATA(pyHsml);

    /* white dwarf data */
    if (!PyDict_Check(wd) ||
        !PyDict_Contains(wd, PyUnicode_FromString("ncells")) ||
        !PyDict_Contains(wd, PyUnicode_FromString("dr"))     ||
        !PyDict_Contains(wd, PyUnicode_FromString("rho")))
    {
        PyErr_SetString(PyExc_ValueError, "wd has do be a dictionary containing at least ncells, dr and rho");
        return nullptr;
    }

    ncells = PyLong_AsLong(PyDict_GetItemString(wd, "ncells"));
    dr = PyFloat_AsDouble(PyDict_GetItemString(wd, "dr"));

    wdRho = (PyArrayObject*)PyDict_GetItemString(wd, "rho");
    if (PyArray_NDIM(wdRho) != 1 || PyArray_DIMS(wdRho)[0] != ncells || PyArray_TYPE(wdRho) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "wd['rho'] has to be of dimensions wd['ncells'] and type double");
        return nullptr;
    }
    rho = (double*)PyArray_DATA(wdRho);

    relaxData(npart, real_pos, mass, hsml, nneighbours, ncells, dr, rho, nsteps);

    pyNewPos = createPyArray2Dim(real_pos, npart, 3);
    free(real_pos);

    return (PyObject*)pyNewPos;
}

static PyMethodDef pysphMethods[] = {
        { "makeTree", makeTree, METH_VARARGS, "" },
        { "relaxData", relaxData, METH_VARARGS, "" },
        { nullptr, nullptr, 0, nullptr }
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pysph", /* m_name */
        nullptr,      /* m_doc */
        -1,                  /* m_size */
        pysphMethods,    /* m_methods */
        nullptr,                /* m_reload */
        nullptr,                /* m_traverse */
        nullptr,                /* m_clear */
        nullptr,                /* m_free */
};

PyMODINIT_FUNC PyInit_pysph(void)
{
    import_array()
    pyTreeType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyTreeType) < 0)
        return nullptr;
    return PyModule_Create(&moduledef);
}


