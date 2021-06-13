#ifndef AREPO_HELPER_LIBS_MAKE_PCOLOR_H
#define AREPO_HELPER_LIBS_MAKE_PCOLOR_H

void check_for_contours(int *neighbours, int *contours, int resolution_x, int resolution_y);
PyObject *make_pcolor(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *make_pcolor_implementation(PyArrayObject *pos, PyArrayObject *quant,
                                     PyArrayObject *axes,
                                     PyArrayObject *boxsizes,
                                     PyArrayObject *resolutions,
                                     PyArrayObject *centers,
                                     int include_neighbours_in_output,
                                     int numthreads);

#endif //AREPO_HELPER_LIBS_MAKE_PCOLOR_H
