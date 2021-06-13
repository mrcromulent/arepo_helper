#include <Python.h>
#include <hdf5.h>
#include <arrayobject.h>
#include <cmath>
#include <vector>

#include "utils.h"

void write_healpix_to_file(PyObject *wdec_3d_dict, const char *fname) {

    const int NTYPES = 6;

    auto pos = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_pos");
    auto ndata_pos = (double *) PyArray_DATA(pos);
    auto u = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_u");
    auto ndata_u = (double *) PyArray_DATA(u);
    auto rho = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_rho");
    auto ndata_rho = (double *) PyArray_DATA(rho);
    auto temp = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_temp");
    auto ndata_temp = (double *) PyArray_DATA(temp);
    auto mass = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_mass");
    auto ndata_mass = (double *) PyArray_DATA(mass);
    auto xnuc = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_xnuc");
    auto ndata_xnuc = (double *) PyArray_DATA(xnuc);
    auto vel = (PyArrayObject *) PyDict_GetItemString(wdec_3d_dict, "ndata_vel");
    auto ndata_vel = (double *) PyArray_DATA(vel);


    const int p = PyArray_DIMS(u)[0];
    int nspecies = PyArray_DIMS(xnuc)[1];

    auto boxsize = PyFloat_AsDouble(PyDict_GetItemString(wdec_3d_dict, "boxsize"));

    // Write to hdf5
//    const char *fname = "bin.dat.ic.hdf5";
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

    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_ThisFile", H5T_NATIVE_INT, buf);
    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_Total", H5T_NATIVE_INT, buf);
    write_vector_to_file(hdf5_headergrp, NTYPES, "NumPart_Total_HighWord", H5T_NATIVE_INT, numpartTotalHighword);
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


    std::vector<const char *> field_names = {"/PartType0/Coordinates", "/PartType0/Velocities",
                                             "/PartType0/NuclearComposition", "/PartType0/InternalEnergy",
                                             "/PartType0/Density", "/PartType0/Temperature", "/PartType0/Masses"};

    std::vector<int> num_columns_v = {3, 3, nspecies, 1, 1, 1, 1};
    std::vector<double *> all_data = {ndata_pos, ndata_vel, ndata_xnuc, ndata_u, ndata_rho, ndata_temp, ndata_mass};
    write_the_hdf5_stuff(hdf5_file, parttype0Group, p, field_names, num_columns_v, all_data);

    auto *particleIDs = static_cast<uint *>(malloc(p * sizeof(uint)));
    for (int l = 0; l < p; l++) {particleIDs[l] = l+1;}

    hsize_t dims[1]; dims[0] = p;
    auto curr_name = "/PartType0/ParticleIDs";
    auto hdf_datatype = H5Tcopy(H5T_NATIVE_UINT);
    auto dataspace_id = H5Screate_simple(1, dims, nullptr);
    auto dataset_id = H5Dcreate1(parttype0Group, curr_name, hdf_datatype, dataspace_id, H5P_DEFAULT);
    write_float_buffer(hdf5_file, p, particleIDs, curr_name, hdf_datatype);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    free(ndata_pos);
    free(ndata_mass);
    free(ndata_rho);
    free(ndata_u);
    free(ndata_vel);
    free(ndata_xnuc);
}

