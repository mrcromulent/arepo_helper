#ifndef AREPO_HELPER_LIBS_MAKE_PCOLOR_H
#define AREPO_HELPER_LIBS_MAKE_PCOLOR_H

void check_for_contours(int *neighbours, int *contours, int resolution_x, int resolution_y);
PyObject *make_pcolor(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *make_pcolor_implementation(PyArrayObject *pos, PyArrayObject *value,
                                     int resolution_x, int resolution_y,
                                     double boxsize_x, double boxsize_y, double boxsize_z,
                                     double center_x, double center_y, double center_z,
                                     int axis0, int axis1,
                                     int nz,
                                     bool include_neighbours_in_output,
                                     int numthreads);

#endif //AREPO_HELPER_LIBS_MAKE_PCOLOR_H
