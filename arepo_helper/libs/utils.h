#ifndef AREPO_HELPER_LIBS_UTILS_H
#define AREPO_HELPER_LIBS_UTILS_H

#include <string>
#include <iostream>
#include <vector>
#include <hdf5.h>


inline long PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {

    // Guards against segfaults
    if(PyArray_API == nullptr) {
        Py_Initialize();
        import_array()
    }

    auto success = PyDict_SetItemString(dict, key, object);
    Py_DECREF(object);
    return success;
}

inline PyArrayObject *createPyArray(double *data, int dim1, int dim2) {
    PyArrayObject *pyData;

    if(PyArray_API == nullptr)
    {
        Py_Initialize();
        import_array()
    }

    if (dim1 == 1 || dim2 == 1) {
        npy_intp dims[1];
        dims[0] = dim1 * dim2;
        pyData = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        memcpy(PyArray_DATA(pyData), data, dim1 * dim2 * sizeof(double));
    } else {
        npy_intp dims[2];
        dims[0] = dim1;
        dims[1] = dim2;
        pyData = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        memcpy(PyArray_DATA(pyData), data, dim1 * dim2 * sizeof(double));
    }

    return pyData;
}
inline PyArrayObject* createPyArray( double *data, int length ) {
    PyArrayObject* pyData;

    npy_intp dims[1];
    dims[0] = length;

    pyData = (PyArrayObject *)PyArray_SimpleNew( 1, dims, NPY_DOUBLE );
    memcpy( PyArray_DATA(pyData), data, length*sizeof(double) );

    return pyData;
}
inline double* resize( double *data, int oldsize, int newsize ) {
    double *newdata;
    int copylength;

    newdata = (double*)malloc( newsize * sizeof(double) );

    copylength = std::min( oldsize, newsize );

    memcpy( newdata, data, copylength * sizeof(double) );
    if (newsize > oldsize) {
        memset( &newdata[oldsize], 0, (newsize-copylength)*sizeof(double) );
    }

    free(data);
    return newdata;
}

inline void print_vector(const std::vector<double> &vec) {
    for (double d : vec)
        std::cout << d << " ";
    std::cout << std::endl;
}

inline double *convert_to_double_star(std::vector<double> vector, int n) {

    auto *retVal = (double *) malloc(n * vector.size() * sizeof(double));

    for (int i = 0; i < vector.size(); i++)
    {
        for (int j = 0; j < n; j++)
        {
            retVal[i * n + j] = vector[i];
        }
    }

    return retVal;
}

static unsigned int Offset = 0;

inline void write_float_buffer_vector(hid_t file_id, unsigned int nof_particles, double *buffer, const char *dataset_name,
                               long hdf_datatype, int num_columns) {
    //identifier
//    herr_t status;

    //open dataset and get dataspace
    hid_t dataset = H5Dopen1(file_id, dataset_name);
    hid_t filespace = H5Dget_space(dataset);

    //velocity file hyperslab
    hsize_t file_offset[2];
    hsize_t file_count[2];
    file_offset[0] = Offset;
    file_offset[1] = 0;
    file_count[0] = nof_particles;
    file_count[1] = num_columns;

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);

    //velocity memory hyperslab
    hsize_t mem[1];
    mem[0] = nof_particles * num_columns;

    hid_t memspace = H5Screate_simple(1, mem, nullptr);
    hsize_t mem_offset[1];
    hsize_t mem_count[1];
    mem_offset[0] = 0;
    mem_count[0] = nof_particles * num_columns;

    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, nullptr, mem_count, nullptr);
    H5Dwrite(dataset, hdf_datatype, memspace, filespace, H5P_DEFAULT, buffer);
    //close handles
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);
}

inline void write_float_buffer(hid_t file_id, unsigned int nof_particles, double *buffer, const char *dataset_name,
                        long hdf_datatype) {
    //identifier
//    herr_t status;

    //open dataset and get dataspace
    hid_t dataset = H5Dopen1(file_id, dataset_name);
    hid_t filespace = H5Dget_space(dataset);

    //file hyperslab
    hsize_t file_offset[1];
    hsize_t file_count[1];
    file_offset[0] = Offset;
    file_count[0] = nof_particles;

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);

    //memory hyperslab
    hsize_t mem[1];
    mem[0] = nof_particles;

    hid_t memspace = H5Screate_simple(1, mem, nullptr);
    hsize_t mem_offset[1];
    hsize_t mem_count[1];
    mem_offset[0] = 0;
    mem_count[0] = nof_particles;

    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, nullptr, mem_count, nullptr);

    //write
    H5Dwrite(dataset, hdf_datatype, memspace, filespace, H5P_DEFAULT, buffer);

    //close handles
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);

}

static void write_float_buffer(hid_t file_id, unsigned int nof_particles, uint *buffer, const char *dataset_name, long hdf_datatype)
{
    //identifier
//    herr_t status;

    //open dataset and get dataspace
    hid_t dataset = H5Dopen1(file_id, dataset_name);
    hid_t filespace = H5Dget_space(dataset);

    //file hyperslab
    hsize_t file_offset[1];
    hsize_t file_count[1];
    file_offset[0] = Offset;
    file_count[0] = nof_particles;

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_offset, nullptr, file_count, nullptr);

    //memory hyperslab
    hsize_t mem[1];
    mem[0] = nof_particles;

    hid_t memspace = H5Screate_simple(1, mem, nullptr);
    hsize_t mem_offset[1];
    hsize_t mem_count[1];
    mem_offset[0] = 0;
    mem_count[0] = nof_particles;

    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, nullptr, mem_count, nullptr);

    //write
    H5Dwrite(dataset, hdf_datatype, memspace, filespace, H5P_DEFAULT, buffer);

    //close handles
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);

}

struct HDFDataField {
    const char* name = "";
    hid_t file_name = {}, group_name = {};
    hid_t datatype = H5T_NATIVE_DOUBLE;
    hid_t dcpl_id = H5P_DEFAULT;
    int p = 0;
    int num_columns = 0;
    double *data = {};

    void write() const {

        long hdf_datatype, dataspace_id, dataset_id;
//        int status;
//        int rank;

        hdf_datatype = H5Tcopy(datatype);

        if (num_columns > 1) { // two dimensional
            hsize_t dims[2]; dims[0] = p; dims[1] = num_columns;
            dataspace_id = H5Screate_simple(2, dims, nullptr);
            dataset_id = H5Dcreate1(group_name, name, hdf_datatype, dataspace_id, dcpl_id);
            write_float_buffer_vector(file_name, p, data, name, hdf_datatype, num_columns);
        } else {
            hsize_t dims[1];
            dims[0] = p;
            dataspace_id = H5Screate_simple(1, dims, nullptr);
            dataset_id = H5Dcreate1(group_name, name, hdf_datatype, dataspace_id, dcpl_id);
            write_float_buffer(file_name, p, data, name, hdf_datatype);
        }

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
};

inline void
write_the_hdf5_stuff(hid_t hdf5_file, hid_t hdf5_group, int num_particles, const std::vector<const char *> &field_names,
                     const std::vector<int> &num_columns_v, const std::vector<double *> &all_data) {

    for(int i = 0; i < num_columns_v.size(); i++) {
        HDFDataField thing;
        thing.name           = field_names[i];
        thing.file_name      = hdf5_file;
        thing.group_name     = hdf5_group;
        thing.p              = num_particles;
        thing.num_columns    = num_columns_v[i];
        thing.data           = all_data[i];
        thing.write();
    }

}

inline void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, double value) {

    auto hdf5_dataspace = H5Screate(H5S_SCALAR);
    auto hdf5_attribute = H5Acreate1(headergrp, name, datatype, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, datatype, &value);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
}

inline void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, int value) {

    auto hdf5_dataspace = H5Screate(H5S_SCALAR);
    auto hdf5_attribute = H5Acreate1(headergrp, name, datatype, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, datatype, &value);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
}

inline void write_vector_to_file(hid_t headergrp, const int ntypes, const char *field_name, hid_t datatype, int *data) {

    hsize_t adim[1] = { static_cast<hsize_t>(ntypes) };

    hid_t hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, nullptr);
    auto hdf5_attribute = H5Acreate1(headergrp, field_name, datatype, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, datatype, data);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

}

#endif //AREPO_HELPER_LIBS_UTILS_H
