#ifndef AREPO_HELPER_LIBS_MAKE_RADIAL_H
#define AREPO_HELPER_LIBS_MAKE_RADIAL_H

#include <Python.h>
#include <armadillo>
#include <arrayobject.h>

PyObject *make_radial(PyObject *self, PyObject *args);
bool inside_cylinder(const arma::vec& A, const arma::vec& B, const arma::vec& P, double radius, arma::vec& Q);
PyArrayObject *make_radial_implementation(PyArrayObject *pos,
                                          PyArrayObject *quant,
                                          PyArrayObject *a,
                                          PyArrayObject *b,
                                          double cylinder_radius,
                                          int nshells);
#endif //AREPO_HELPER_LIBS_MAKE_RADIAL_H
