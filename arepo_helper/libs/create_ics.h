#ifndef AREPO_HELPER_LIBS_CREATE_ICS_H
#define AREPO_HELPER_LIBS_CREATE_ICS_H

#include <Python.h>
#include <arrayobject.h>

typedef struct {
    double x, y, z, r;
} t_point;

int compare_points(const void *o1, const void *o2);
PyObject *create_particles_cube(PyObject *self, PyObject *args);
PyObject *create_particles_healpix(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *create_particles_fill_grid(PyObject *self, PyObject *args, PyObject *kwargs);

#endif //AREPO_HELPER_LIBS_CREATE_ICS_H
