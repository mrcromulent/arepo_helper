#include <Python.h>
#include <armadillo>
#include <arrayobject.h>

#include "make_radial.h"
#include "make_wdec.h"


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

PyObject * make_radial_check_inputs(PyObject *args,
                                    PyArrayObject *pos, PyArrayObject *data,
                                    int nshells) {

    if (!PyArg_ParseTuple(args, "O!O!i:make_radial(pos, data, nshells)",
                                     &PyArray_Type, &pos, &PyArray_Type, &data, &nshells)) {
        return nullptr;
    }

    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    if (PyArray_NDIM(data) != 1 || PyArray_TYPE(data) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "data has to be of dimension [n] and type double");
        return nullptr;
    }
    return args;
}

//PyObject *make_radial(PyArrayObject *pos, PyArrayObject *data) {
PyObject *make_radial(PyObject *self, PyObject *args) {

    PyArrayObject *pos, *data;
    int nshells = 0;

    if (!PyArg_ParseTuple(args, "O!O!i:make_radial(pos, data, nshells)",
                          &PyArray_Type, &pos, &PyArray_Type, &data, &nshells)) {
        return nullptr;
    }

    if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [n,3] and type double");
        return nullptr;
    }

    if (PyArray_NDIM(data) != 1 || PyArray_TYPE(data) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "data has to be of dimension [n] and type double");
        return nullptr;
    }

    auto pyProfile = make_radial_implementation(pos, data, nshells);

    auto value_data = (double *) PyArray_DATA(pyProfile);
    for (int i = 0; i < nshells; i++) {
        std::cout << value_data[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < nshells; i++) {
        std::cout << value_data[nshells + i] << " ";
    }
    std::cout << std::endl;

    return PyArray_Return(pyProfile);
}

PyArrayObject *make_radial_implementation(PyArrayObject *pos, PyArrayObject *data, int nshells) {
    int npart = PyArray_DIMS(pos)[0];

    // Guards against segfaults
    if(PyArray_API == nullptr) {
        Py_Initialize();
        import_array()
    }

    npy_intp dims[2] = {2, nshells};
    auto pyProfile = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    std::cout << "made py profile" << "\n";
    auto data_pos = (double *) PyArray_DATA(pos);
    std::cout << "got data_pos" << "\n";
    auto data_data = (double *) PyArray_DATA(data);
    std::cout << "got data_data" << "\n";

    std::cout << "got npart" << "\n";

    // Start the clock
    clock_t start = clock();

    std::cout << "clock started" << "\n";

    auto count = (int *) malloc(nshells * sizeof(int));
    auto profile = (double *) PyArray_DATA(pyProfile);
    memset(profile, 0, 2 * nshells * sizeof(double));
    memset(count, 0, nshells * sizeof(int));

    double boxsize = 1e10;
    double cylinder_r = 0.01 * boxsize;
    arma::vec A = {boxsize/2, boxsize/2, boxsize/2};
    arma::vec B = {boxsize, boxsize/2, boxsize/2};
    double d_total = arma::norm(B - A);

    arma::vec Q;
    arma::vec P;
    double px, py, pz, d_i;
    bool inside;

    std::cout << "d_total = " << d_total << std::endl;

    data_pos = (double *) PyArray_DATA(pos);
    data_data = (double *) PyArray_DATA(data);


    for (int part = 0; part < npart; part++) {
        px = *data_pos;
        data_pos = (double *) ((char *) data_pos + PyArray_STRIDES(pos)[1]);
        py = *data_pos;
        data_pos = (double *) ((char *) data_pos + PyArray_STRIDES(pos)[1]);
        pz = *data_pos;
        data_pos = (double *) ((char *) data_pos - 2 * PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);

        auto d = *data_data;
        data_data = (double *) ((char *) data_data + PyArray_STRIDES(data)[0]);

        P = {px, py, pz};

        inside = inside_cylinder(A, B, P, cylinder_r, Q);

        if (inside) {
//            std::cout << "Just about to perform subtraction" << "\n";
            d_i = arma::norm(Q - A);

            int shell = floor(d_i / d_total * nshells);
//            std::cout << shell << std::endl;
            profile[shell]  += d;
            count[shell]    += 1;
        }
    }

    for (int shell = 0; shell < nshells; shell++) {
        profile[nshells + shell] = (double) shell / (double) nshells * d_total;

        if (count[shell] > 0) {
            profile[shell] = profile[shell] / count[shell];
        }
    }



    free(count);
    printf("Calculation took %gs\n", ((double) clock() - (double) start) / CLOCKS_PER_SEC);

    return pyProfile;

}

// Python Module definition
static PyMethodDef radial_methods[] = {
        {"make_radial", make_radial, METH_VARARGS, ""},
        {nullptr, nullptr, 0,nullptr}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "radial_pierre",
        nullptr,
        -1,
        radial_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_radial_pierre(void) {
    import_array();

    return PyModule_Create(&moduledef);
}