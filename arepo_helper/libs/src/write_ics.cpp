#include <Python.h>
#include <hdf5.h>
#include <arrayobject.h>
#include <cmath>
#include <vector>

#include "../headers/utils.h"
#include "../headers/const.h"

void write_dict_to_hdf5(PyObject *dict, const char *fname, double boxsize) {

    const int NTYPES = 6;

    auto pos = (PyArrayObject *) PyDict_GetItemString(dict, f[N::COORDINATES]);
    auto ndata_pos = (double *) PyArray_DATA(pos);
    auto u = (PyArrayObject *) PyDict_GetItemString(dict, f[N::INTERNALENERGY]);
    auto ndata_u = (double *) PyArray_DATA(u);
    auto rho = (PyArrayObject *) PyDict_GetItemString(dict, f[N::DENSITY]);
    auto ndata_rho = (double *) PyArray_DATA(rho);
    auto mass = (PyArrayObject *) PyDict_GetItemString(dict, f[N::MASSES]);
    auto ndata_mass = (double *) PyArray_DATA(mass);
    auto pressure = (PyArrayObject *) PyDict_GetItemString(dict, f[N::PRESSURE]);
    auto ndata_pres = (double *) PyArray_DATA(pressure);
    auto xnuc = (PyArrayObject *) PyDict_GetItemString(dict, f[N::NUCLEARCOMPOSITION]);
    auto ndata_xnuc = (double *) PyArray_DATA(xnuc);
    auto vel = (PyArrayObject *) PyDict_GetItemString(dict, f[N::VELOCITIES]);
    auto ndata_vel = (double *) PyArray_DATA(vel);


    const int p = PyArray_DIMS(u)[0];
    int nspecies = PyArray_DIMS(xnuc)[1];


    // Write to hdf5
    auto hdf5_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    auto hdf5_headergrp = H5Gcreate1(hdf5_file, "/Header", 0);
    auto parttype0Group = H5Gcreate1(hdf5_file, "/PartType0", 0);


    int buf[NTYPES]                     = {p, 0, 0, 0, 0, 0};
    int numpartTotalHighword[NTYPES]    = {0, 0, 0, 0, 0, 0};
    double massTable[NTYPES]            = {0, 0, 0, 0, 0, 0};
    hsize_t adim[1]                     = {NTYPES};

    hid_t hdf5_dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(hdf5_dataspace, 1, adim, nullptr);
    hid_t hdf5_attribute = H5Acreate1(hdf5_headergrp, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, massTable);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);

    // Write vector quantities
    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_ThisFile", H5T_NATIVE_INT, buf);
    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_Total", H5T_NATIVE_INT, buf);
    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_Total_HighWord", H5T_NATIVE_INT, numpartTotalHighword);

    // Write scalar quantities
    write_scalar_to_file(hdf5_headergrp, "Time", H5T_NATIVE_DOUBLE, 0.0);
    write_scalar_to_file(hdf5_headergrp, "BoxSize", H5T_NATIVE_DOUBLE, boxsize);
    write_scalar_to_file(hdf5_headergrp, "Redshift", H5T_NATIVE_DOUBLE, 0.0);
    write_scalar_to_file(hdf5_headergrp, "Omega0", H5T_NATIVE_DOUBLE, 0.0);
    write_scalar_to_file(hdf5_headergrp, "OmegaLambda", H5T_NATIVE_DOUBLE, 0.0);
    write_scalar_to_file(hdf5_headergrp, "HubbleParam", H5T_NATIVE_DOUBLE, 0.0);
    write_scalar_to_file(hdf5_headergrp, "NumFilesPerSnapshot", H5T_NATIVE_INT, 1);
    write_scalar_to_file(hdf5_headergrp, "Flag_DoublePrecision", H5T_NATIVE_INT, 1);
    write_scalar_to_file(hdf5_headergrp, "Flag_Feedback", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_Sfr", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_Cooling", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_StellarAge", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_Metals", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_IC_info", H5T_NATIVE_INT, 0);
    write_scalar_to_file(hdf5_headergrp, "Flag_Entropy_ICs", H5T_NATIVE_INT, 0);


    std::vector<const char *> field_names = {"/PartType0/Coordinates",
                                             "/PartType0/Velocities",
                                             "/PartType0/NuclearComposition",
                                             "/PartType0/InternalEnergy",
                                             "/PartType0/Density",
                                             "/PartType0/Masses",
                                             "/PartType0/Pressure"};

    std::vector<int> num_columns_v = {3, 3, nspecies, 1, 1, 1, 1};
    std::vector<double *> all_data = {ndata_pos, ndata_vel, ndata_xnuc, ndata_u, ndata_rho, ndata_mass, ndata_pres};
    write_the_hdf5_stuff(hdf5_file, parttype0Group, p, field_names, num_columns_v, all_data);

    free(ndata_pos);
    free(ndata_mass);
    free(ndata_rho);
    free(ndata_u);
    free(ndata_vel);
    free(ndata_xnuc);
    free(ndata_pres);
}

