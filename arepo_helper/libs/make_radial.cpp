#include <Python.h>
#include <armadillo>
#include <arrayobject.h>

#include "make_radial.h"


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

PyObject *make_radial(PyObject *self, PyObject *args) {

    PyArrayObject *pos, *quant, *a, *b;
    double cylinder_radius;
    int nshells = 0;

    if (!PyArg_ParseTuple(args, "O!O!O!O!di:make_radial(pos, quant, a, b, cylinder_radius, nshells)",
                          &PyArray_Type, &pos, &PyArray_Type, &quant,
                          &PyArray_Type, &a, &PyArray_Type, &b,
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

PyArrayObject *make_radial_implementation(PyArrayObject *pos,
                                          PyArrayObject *quant,
                                          PyArrayObject *a,
                                          PyArrayObject *b,
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
    double a_x  = * (double *) ((char *) PyArray_DATA(a) + 0 * PyArray_STRIDES(a)[0]);
    double a_y  = * (double *) ((char *) PyArray_DATA(a) + 1 * PyArray_STRIDES(a)[0]);
    double a_z  = * (double *) ((char *) PyArray_DATA(a) + 2 * PyArray_STRIDES(a)[0]);
    double b_x  = * (double *) ((char *) PyArray_DATA(b) + 0 * PyArray_STRIDES(a)[0]);
    double b_y  = * (double *) ((char *) PyArray_DATA(b) + 1 * PyArray_STRIDES(a)[0]);
    double b_z  = * (double *) ((char *) PyArray_DATA(b) + 2 * PyArray_STRIDES(a)[0]);

    auto data_pos = (double *) PyArray_DATA(pos);
    auto data_data = (double *) PyArray_DATA(quant);

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
    double px, py, pz, d_i;
    bool inside;

    data_pos = (double *) PyArray_DATA(pos);
    data_data = (double *) PyArray_DATA(quant);


    for (int part = 0; part < npart; part++) {
        px = *data_pos;
        data_pos = (double *) ((char *) data_pos + PyArray_STRIDES(pos)[1]);
        py = *data_pos;
        data_pos = (double *) ((char *) data_pos + PyArray_STRIDES(pos)[1]);
        pz = *data_pos;
        data_pos = (double *) ((char *) data_pos - 2 * PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);

        auto d = *data_data;
        data_data = (double *) ((char *) data_data + PyArray_STRIDES(quant)[0]);

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
        "arepo_radial",
        nullptr,
        -1,
        radial_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_arepo_radial(void) {
    import_array();

    return PyModule_Create(&moduledef);
}