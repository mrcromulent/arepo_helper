#ifndef AREPO_HELPER_LIBS_CREATE_ICS_H
#define AREPO_HELPER_LIBS_CREATE_ICS_H

#include <Python.h>
#include <arrayobject.h>

typedef struct {
    double x, y, z, r;
} t_point;

int compare_points(const void *o1, const void *o2);
PyObject *create_particles_cube(PyObject *self, PyObject *args);
PyObject *create_particles_fill_grid(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject *convert_to_healpix_implementation(PyObject *wdec_dict,
                                            int nspecies,
                                            double boxsize,
                                            PyArrayObject *centers_py,
                                            int makebox,
                                            int randomizeshells,
                                            int randomizeradii,
                                            double pmass);

#endif //AREPO_HELPER_LIBS_CREATE_ICS_H
