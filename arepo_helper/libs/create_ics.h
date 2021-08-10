#ifndef AREPO_HELPER_LIBS_CREATE_ICS_H
#define AREPO_HELPER_LIBS_CREATE_ICS_H

#include <Python.h>
#include <arrayobject.h>


PyObject *convert_to_healpix_implementation(PyObject *dict_in, double boxsize, PyObject *centers_py,
                                            int randomizeshells, int randomizeradii, double pmass);


PyArrayObject *create_particles_fill_grid_implementation(PyArrayObject* pyPos,
                                                    double boxsize,
                                                    int boxres);

PyObject *add_grid_particles_implementation(PyObject *dict,
                                            double boxsize,
                                            int boxres,
                                            double grid_pres,
                                            double grid_density,
                                            PyObject *grid_xnuc);

PyArrayObject *create_particles_fill_grid(PyObject *self, PyObject *args);

PyObject *convert_to_healpix(PyObject *self, PyObject *args, PyObject *kwargs);

PyObject *add_grid_particles(PyObject *self, PyObject *args, PyObject *kwargs);

// Python Module Definition
static PyMethodDef create_ics_methods[] = {
        {"convert_to_healpix", (PyCFunction) convert_to_healpix,METH_VARARGS|METH_KEYWORDS,""},
        {"create_particles_fill_grid", (PyCFunction) create_particles_fill_grid, METH_VARARGS,""},
        {"add_grid_particles", (PyCFunction) add_grid_particles, METH_VARARGS|METH_KEYWORDS,""},
        {nullptr,nullptr,0,nullptr}
};

static struct PyModuleDef moduledef_create_ics = {
        PyModuleDef_HEAD_INIT,
        "create_ics",
        nullptr,
        -1,
        create_ics_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};


#endif //AREPO_HELPER_LIBS_CREATE_ICS_H
