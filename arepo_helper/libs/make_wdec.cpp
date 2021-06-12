#include <Python.h>
#include <hdf5.h>
#include <arrayobject.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <vector>


#include "mersenne.h"
#include "pix2vec_ring.h"
#include "spline.h"
#include "make_wdec.h"
#include "const.h"
#include "utils.h"




PyObject *make_wdec_new(const char *wdec_dir,
                   double boxsize,
                   int nspecies,
                   bool makebox,
                   bool randomizeshells,
                   bool randomizeradii,
                   double pmass) {

    seed();


    const int NTYPES = 6;
    double gamma = 5.0 / 3.0;
    double num_base_pixels = 12;

    double grid_energy  = 0;
    double grid_density = 1e-4;
    double grid_temp    = 1500;
    int boxres          = 32;
    int npart_max       = 1e8;


    auto wd_results = WdecResults(wdec_dir);
    int npoints         = wd_results.n_points;
    auto r_extended     = wd_results.radii;
    auto mr_extended    = wd_results.mr;
    auto temp_extended  = wd_results.temp;
    auto rho_extended   = wd_results.density;
    auto pres_extended  = wd_results.pres;
    auto xo_extended    = wd_results.xo;
    auto xc_extended    = wd_results.xc;



    //
    std::vector<double> e_extended;
    e_extended.reserve(rho_extended.size());
    for (int i = 0; i < rho_extended.size(); i++) {
        e_extended.push_back(pres_extended[i] / (rho_extended[i] * (gamma - 1)));
    }


    double mtot = mr_extended[mr_extended.size() - 1];
    printf("Total mass: %g solar masses\n", mtot / msol);

    // convert to double *
    auto data_r     = convert_to_double_star(r_extended, 1);
    auto data_mr    = convert_to_double_star(mr_extended, 1);



    tk::spline spline_carbon(r_extended, xc_extended);
    tk::spline spline_oxygen(r_extended, xo_extended);
    tk::spline spline_u(r_extended, e_extended);
    tk::spline spline_rho(r_extended, rho_extended);
    tk::spline spline_temp(r_extended, temp_extended);

    auto three_vector_data_size     = 3 * npart_max * sizeof(double);
    auto scalar_data_size           = npart_max * sizeof(double);
    auto ndata_pos                  = (double *) malloc(three_vector_data_size);
    auto ndata_vel                  = (double *) malloc(three_vector_data_size);
    auto ndata_mass                 = (double *) malloc(scalar_data_size);
    auto ndata_temp                 = (double *) malloc(scalar_data_size);
    auto ndata_rho                  = (double *) malloc(scalar_data_size);
    auto ndata_u                    = (double *) malloc(scalar_data_size);
    auto ndata_xnuc                 = (double *) malloc(nspecies * npart_max * sizeof(double));

    int max_nside = 15;
    double inner_radius, outer_radius, shell_radius;
    double interior_mass;
    double width = 0;
    double phis[3], hp_vector[3];
    int n_pix;
    long nside;
    double n2;

    int p = 0;

    for (int i = 0; i < (npoints - 1); i++) {


        inner_radius = data_r[i];
        outer_radius = data_r[i + 1];
        shell_radius = 0.5 * (inner_radius + outer_radius);
        interior_mass = data_mr[i];

        // TODO pmass might not be the best way to achieve this
        n2      = sqrt(interior_mass / pmass / num_base_pixels); /* increases with index */
        nside   = floor(n2 + 0.5);
        if (nside >= max_nside) { nside = max_nside;}
        if (nside < 1) {nside = 1;}
        n_pix = static_cast<int>(num_base_pixels * pow(nside, 2));


        if (randomizeshells) {
            for (double &phi : phis) { phi = randMT() * 2. * M_PI; }
        }

        for (int ipring = 0; ipring < n_pix; ipring++) {
            // Renders cartesian vector coordinates of the normal pixel center given the pixel number ipring and
            // a map resolution parameter nside
            pix2vec_ring(nside, ipring, hp_vector);

            if (randomizeradii)
                shell_radius += 0.1 * width * (randMT() - 0.5);

            if (randomizeshells) {

                double phi1 = phis[0], phi2 = phis[1], phi3 = phis[2];

                double x = hp_vector[0];
                double y = cos(phi1) * hp_vector[1] - sin(phi1) * hp_vector[2];
                double z = sin(phi1) * hp_vector[1] + cos(phi1) * hp_vector[2];

                double x2 = cos(phi2) * x + sin(phi2) * z;
                double y2 = y;
                double z2 = -sin(phi2) * x + cos(phi2) * z;

                // TODO: generalise the centering policy
                ndata_pos[p * 3 + 0] = shell_radius * (cos(phi3) * x2 - sin(phi3) * y2) + boxsize / 2;
                ndata_pos[p * 3 + 1] = shell_radius * (sin(phi3) * x2 + cos(phi3) * y2) + boxsize / 2;
                ndata_pos[p * 3 + 2] = shell_radius * z2 + boxsize / 2;
            } else {
                for (int k = 0; k < 3; k++)
                    ndata_pos[p * 3 + k] = hp_vector[k] * shell_radius + boxsize / 2;
            }


            // blah
            ndata_mass[p]   = pmass;
            ndata_rho[p]    = spline_rho(shell_radius);
            ndata_u[p]      = spline_u(shell_radius);
            ndata_temp[p]   = spline_temp(shell_radius);
            double xCHere   = spline_carbon(shell_radius);
            double xOHere   = spline_oxygen(shell_radius);

            for (int k = 0; k < nspecies; k++) {
                ndata_xnuc[p * nspecies + k] = 0;
            }

            if (xCHere + xOHere >= 1) {xOHere = 1 - xCHere;}
            ndata_xnuc[p * nspecies + 0] = 1 - xCHere - xOHere; // He
            ndata_xnuc[p * nspecies + 1] = xCHere; // C
            ndata_xnuc[p * nspecies + 2] = xOHere; // O

            // Set a density floor
            if (ndata_rho[p] < grid_density) {
                ndata_rho[p] = grid_density;
            }

            // Iterate the particle counter
            p++;
        }

    }

    memset(ndata_vel, 0, 3 * p * sizeof(double));
    auto *particleIDs = static_cast<uint *>(malloc(p * sizeof(uint)));
    for (int l = 0; l < p; l++) {particleIDs[l] = l+1;}

    printf("Created %d particles.\n", p);


    if (makebox) {
        // Nuclear composition of grid particles
        auto grid_xnuc = (double *) malloc(nspecies * sizeof(double));
        for (int k = 0; k < nspecies; k++) {grid_xnuc[k] = 0;}
        grid_xnuc[0] = 1;

        int nbox            = pow(boxres, 3);
        double boxcellsize  = boxsize / (double) boxres;
        double volboxcell   = pow(boxcellsize, 3);


        char *box = static_cast<char *>(malloc(nbox));
        memset(box, 0, nbox);
        for (int i = 0; i < p; i++) {
            int ix = floor(ndata_pos[i * 3 + 0] / boxcellsize);
            int iy = floor(ndata_pos[i * 3 + 1] / boxcellsize);
            int iz = floor(ndata_pos[i * 3 + 2] / boxcellsize);

            box[iz * boxres * boxres + iy * boxres + ix] = 1;
        }

        int boxCount = 0;
        for (int i = 0; i < nbox; i++) {
            int ix = i % boxres;
            int iy = (i / boxres) % boxres;
            int iz = i / (boxres * boxres);

            double x = (ix + 0.5) * boxcellsize;
            double y = (iy + 0.5) * boxcellsize;
            double z = (iz + 0.5) * boxcellsize;

            if (!box[i]) {
                ndata_pos[p * 3 + 0]    = x;
                ndata_pos[p * 3 + 1]    = y;
                ndata_pos[p * 3 + 2]    = z;
                ndata_rho[p]            = grid_density;
                ndata_mass[p]           = grid_density * volboxcell;
                ndata_temp[p]           = grid_temp;
                ndata_u[p]              = grid_energy;

                for (int k = 0; k < nspecies; k++) {
                    ndata_xnuc[p * nspecies + k] = grid_xnuc[k];
                }

                p++;
                boxCount++;
            }
        }
        free(box);
        free(grid_xnuc);
        printf("Added %d box cells.\n", boxCount);
    }

    // Free all unused memory space
    ndata_xnuc  = (double *) realloc(ndata_xnuc, nspecies * p * sizeof(double));
    ndata_pos   = (double *) realloc(ndata_pos, 3 * p * sizeof(double));
    ndata_vel   = (double *) realloc(ndata_vel, 3 * p * sizeof(double));
    ndata_temp  = (double *) realloc(ndata_temp, p * sizeof(double));
    ndata_mass  = (double *) realloc(ndata_mass, p * sizeof(double));
    ndata_rho   = (double *) realloc(ndata_rho, p * sizeof(double));
    ndata_u     = (double *) realloc(ndata_u, p * sizeof(double));

    // Write to hdf5
    const char *fname = "bin.dat.ic.hdf5";
    auto hdf5_file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    auto hdf5_headergrp = H5Gcreate1(hdf5_file, "/Header", 0);
    auto parttype0Group = H5Gcreate1(hdf5_file, "/PartType0", 0);


    int buf[NTYPES]                     = {p, 0, 0, 0, 0, 0};
    int numpartTotalHighword[NTYPES]    = {0, 0, 0, 0, 0, 0};
    double massTable[NTYPES]            = {0, 0, 0, 0, 0, 0};
    hsize_t adim[1]                     = { NTYPES };

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

    hsize_t dims[1]; dims[0] = p;
    auto curr_name = "/PartType0/ParticleIDs";
    auto hdf_datatype = H5Tcopy(H5T_NATIVE_UINT);
    auto dataspace_id = H5Screate_simple(1, dims, nullptr);
    auto dataset_id = H5Dcreate1(parttype0Group, curr_name, hdf_datatype, dataspace_id, H5P_DEFAULT);
    write_float_buffer(hdf5_file, p, particleIDs, curr_name, hdf_datatype);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    auto dict = PyDict_New();
    PyDict_SetStolenItem(dict, "pos", (PyObject *) createPyArray(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, "mass", (PyObject *) createPyArray(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, "rho", (PyObject *) createPyArray(ndata_rho, p, 1));
    free(ndata_rho);
    PyDict_SetStolenItem(dict, "u", (PyObject *) createPyArray(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, "temp", (PyObject *) createPyArray(ndata_temp, p, 1));
    free(ndata_temp);
    PyDict_SetStolenItem(dict, "vel", (PyObject *) createPyArray(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, "xnuc", (PyObject *) createPyArray(ndata_xnuc, p, nspecies));
    free(ndata_xnuc);
    PyDict_SetStolenItem(dict, "count", (PyObject *) PyLong_FromLong(p));

    return dict;
}
