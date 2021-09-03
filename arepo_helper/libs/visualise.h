#ifndef AREPO_HELPER_LIBS_VISUALISE_H
#define AREPO_HELPER_LIBS_VISUALISE_H

#include <Python.h>
#include <armadillo>
#include <arrayobject.h>

void check_for_contours(const int *neighbours, int *contours, int resolution_x, int resolution_y);

bool inside_cylinder(const arma::vec& A, const arma::vec& B, const arma::vec& P, double radius, arma::vec& Q);

PyArrayObject *make_radial_implementation(PyArrayObject *pos,
                                          PyArrayObject *quant,
                                          PyObject *a,
                                          PyObject *b,
                                          double cylinder_radius,
                                          int nshells);

PyObject *make_pcolor_implementation(PyArrayObject *pos, PyArrayObject *quant,
                                     PyObject *axes,
                                     PyObject *boxsizes,
                                     PyObject *resolutions,
                                     PyObject *centers,
                                     int include_neighbours_in_output,
                                     int numthreads);

PyObject *make_radial(PyObject *self, PyObject *args);

PyObject *make_pcolor(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject *get_indices_of_neighbours(PyObject *self, PyObject *args);

// Python Module definition
static PyMethodDef vis_methods[] = {
        {"make_radial", make_radial, METH_VARARGS, ""},
        {"make_pcolor", (PyCFunction) make_pcolor, METH_VARARGS | METH_KEYWORDS, ""},
        {"get_indices_of_neighbours", get_indices_of_neighbours, METH_VARARGS, ""},
        {nullptr, nullptr, 0,nullptr}
};

static struct PyModuleDef moduledef_arepo_vis = {
        PyModuleDef_HEAD_INIT,
        "arepo_vis",
        nullptr,
        -1,
        vis_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

#endif //AREPO_HELPER_LIBS_VISUALISE_H
