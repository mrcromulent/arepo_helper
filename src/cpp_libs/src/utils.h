#ifndef AREPO_HELPER_LIBS_UTILS_H
#define AREPO_HELPER_LIBS_UTILS_H

#include <Python.h>
#include <arrayobject.h>
#include <vector>
#include <hdf5.h>

double mr_diff(double q1, double q2);
std::vector<double> convert_to_std_vector(PyArrayObject *obj);
bool contains(std::vector<const char*>& v, const char *key);
double* resize(double *data, int old_size, int new_size, double fill_val = 0.0);
long PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object);
void print(PyObject * to_print);
PyArrayObject *create_numpy_array(double *data, int dim1, int dim2);
PyArrayObject *create_numpy_array(double *data, int length);
PyArrayObject *create_numpy_array(int *data, int length);
double *convert_to_double_star(std::vector<double> vector, int n);
void write_float_buffer_vector(hid_t file_id, unsigned int nof_particles, double *buffer, const char *dataset_name,
                               long hdf_datatype, int num_columns);
void write_float_buffer(hid_t file_id, unsigned int nof_particles, double *buffer, const char *dataset_name,
                        long hdf_datatype);
void write_float_buffer(hid_t file_id, unsigned int nof_particles, uint *buffer, const char *dataset_name, long hdf_datatype);
void write_the_hdf5_stuff(hid_t hdf5_file, hid_t hdf5_group,
                          int num_particles,
                          const std::vector<const char *> &field_names,
                          const std::vector<int> &num_columns_v,
                          const std::vector<double *> &all_data);

void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, double value);
void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, int value);
void write_vector_to_file(hid_t headergrp, const int ntypes, const char *field_name, hid_t datatype, int *data);


#endif //AREPO_HELPER_LIBS_UTILS_H
