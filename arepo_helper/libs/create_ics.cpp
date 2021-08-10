#include <Python.h>
#include <arrayobject.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "mersenne.h"
#include "pix2vec_ring.h"
#include "utils.h"
#include "create_ics.h"
#include "spline.h"
#include "const.h"


PyArrayObject *create_particles_fill_grid(PyObject *self, PyObject *args) {

    PyArrayObject *pyPos;
    double boxsize;
    int boxres;

    if (!PyArg_ParseTuple(args, "O!di:create_particles_fill_grid(pos, boxsize, boxres)",
                          &PyArray_Type, &pyPos, &boxsize, &boxres)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "Coordinates has to be of dimensions [N,3] and type double");
        return nullptr;
    }

    auto pos_array = create_particles_fill_grid_implementation(pyPos, boxsize, boxres);

    return pos_array;
}

PyObject *add_grid_particles_implementation(PyObject *dict,
                                            double boxsize,
                                            int boxres,
                                            double grid_pres,
                                            double grid_density,
                                            PyObject *grid_xnuc) {

    double gamma        = 1.4;
    double boxcellsize  = boxsize / (double) boxres;
    double volboxcell   = pow(boxcellsize, 3);
    int nspecies        = PyList_Size(grid_xnuc);


    auto curr_pos       = (PyArrayObject *) PyDict_GetItemString(dict, f[N::COORDINATES]);


    // Find the positions of box particles
    auto box_part_pos   = create_particles_fill_grid_implementation(curr_pos, boxsize, boxres);
    int p_added         = PyArray_DIMS(box_part_pos)[0];
    int curr_particles  = PyArray_DIMS(curr_pos)[0];
    int p               = curr_particles + p_added;

    std::vector<const char *> fields = {f[N::DENSITY], f[N::MASSES], f[N::INTERNALENERGY], f[N::PRESSURE]};
    std::vector<const char *> exclusions = {f[N::COORDINATES], f[N::NUCLEARCOMPOSITION]};
    std::vector<double> vals = {grid_density,
                                grid_density * volboxcell,
                                grid_pres / (grid_density * (gamma - 1)),
                                grid_pres};

    PyObject *key, *keys;
    keys = PyDict_Keys(dict);
    auto ret_dict = PyDict_New();

    for (size_t i = 0; i < PyList_Size(keys); i++) {

        key         = PyList_GetItem(keys, i);
        const char* str = PyUnicode_AsUTF8(key);


        auto item   = (PyArrayObject *) PyDict_GetItem(dict, key);
        auto data   = (double *) PyArray_DATA(item);
        int new_size, old_size;


        if (PyArray_NDIM(item) == 2) {
            new_size = p * (int) PyArray_DIMS(item)[1];
            old_size = curr_particles * (int)  PyArray_DIMS(item)[1];
        } else {
            new_size = p;
            old_size = curr_particles;
        }

        if (contains(fields, str)) {
            int num = 0;
            for (const char *stri : fields) {
                if(strcmp(stri, str) == 0) {
                    data = resize(data, old_size, new_size, vals[num]);
                }
                num++;
            }

        } else { // Fill with zeros if no specific fill value given
            data = resize(data, old_size, new_size);
        }

        if (PyArray_NDIM(item) == 2) {
//            PyDict_SetItem(ret_dict, key, (PyObject *) create_numpy_array(data, p, PyArray_DIMS(item)[1]));
            PyDict_SetStolenItem(ret_dict, str, (PyObject *) create_numpy_array(data, p, PyArray_DIMS(item)[1]));
        } else {
//            PyDict_SetItem(ret_dict, key, (PyObject *) create_numpy_array(data, p));
            PyDict_SetStolenItem(ret_dict, str, (PyObject *) create_numpy_array(data, p));
        }

        Py_DECREF(key);
    }

    Py_DECREF(keys);


    int idx;
    auto data_pos_py    = (PyArrayObject *) PyDict_GetItemString(dict, f[N::COORDINATES]);
    auto data_xnuc_py   = (PyArrayObject *) PyDict_GetItemString(dict, f[N::NUCLEARCOMPOSITION]);
    auto s_pos          = PyArray_STRIDES(data_pos_py);
    auto s_xnuc         = PyArray_STRIDES(data_xnuc_py);
    auto ndata_pos      = (double *) malloc(3 * p * sizeof(double));
    auto ndata_xnuc     = (double *) malloc(nspecies * p * sizeof(double));

    for (int i = 0; i < curr_particles; i++) {
        // Fill in the grid xnucs
        for (int j = 0; j < nspecies; j++) {
            ndata_xnuc[i * nspecies + j] = * (double *) ((char *) PyArray_DATA(data_xnuc_py) + i * s_xnuc[0] + j * s_xnuc[1]);
        }

        // Copy in the new positions
        for (int j = 0; j < 3; j++) {
            ndata_pos[i * 3 + j] = * (double *) ((char *) PyArray_DATA(data_pos_py) + i * s_pos[0] + j * s_pos[1]);
        }
    }


    for (int i = 0; i < p_added; i++) {
        idx = curr_particles + i;

        // Fill in the grid xnucs
        for (int j = 0; j < nspecies; j++) {
            ndata_xnuc[idx * nspecies + j] = PyFloat_AsDouble(PyList_GetItem(grid_xnuc, j));
        }

        // Copy in the new positions
        for (int j = 0; j < 3; j++) {
            ndata_pos[idx * 3 + j] = * (double *) ((char *) PyArray_DATA(box_part_pos) + i * s_pos[0] + j * s_pos[1]);
        }
    }

    PyDict_SetStolenItem(ret_dict, f[N::COORDINATES], (PyObject *) create_numpy_array(ndata_pos, p, 3));
    PyDict_SetStolenItem(ret_dict, f[N::NUCLEARCOMPOSITION], (PyObject *) create_numpy_array(ndata_xnuc, p, nspecies));

    free(ndata_pos);
    free(ndata_xnuc);

    return ret_dict;
}

PyObject *add_grid_particles(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyObject *dict;
    double boxsize;
    int boxres          = 32;
    double grid_pres    = 2e6;
    double grid_density = 1e-4;
    PyObject *grid_xnuc = nullptr;

    const char *kwlist[] = {"dict", "boxsize", "boxres", "grid_pres", "grid_density", "grid_xnuc", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Od|iddO:add_grid_particles(dict, boxsize, [boxres, grid_pres, grid_density, grid_xnuc])",
                                keywords, &dict, &boxsize, &boxres, &grid_pres, &grid_density, &grid_xnuc)) {
        return nullptr;
    }

    // Create Helium background if no grid_xnuc given
    if (grid_xnuc == nullptr) {
        grid_xnuc = Py_BuildValue("[ddddd]", 1.0, 0.0, 0.0, 0.0, 0.0);
    }

    auto expanded_dict = add_grid_particles_implementation(dict, boxsize, boxres, grid_pres, grid_density, grid_xnuc);

    return expanded_dict;
}

PyArrayObject *create_particles_fill_grid_implementation(PyArrayObject* pyPos,
                                                    double boxsize,
                                                    int boxres) {
    auto pos = (double *) PyArray_DATA(pyPos);
    int npart = PyArray_DIMS(pyPos)[0];

    printf("%d particles found.\n", npart);

    double cellsize = boxsize / boxres;
    double boxx = boxsize;
    double boxy = boxsize;
    double boxz = boxsize;

    int nx = floor(boxx / cellsize);
    int ny = floor(boxy / cellsize);
    int nz = floor(boxz / cellsize);
    int ncells = nx * ny * nz;

    printf("Building grid with %d x %d x %d cells.\n", nx, ny, nz);


    auto s_pos = PyArray_STRIDES(pyPos);
    int grid_index;
    double px, py, pz;
    int ix, iy, iz;

    // Initialise a grid which holds the number of particles within each cell
    // Zero initially
    auto grid = static_cast<int *>(malloc(ncells * sizeof(int)));
    memset(grid, 0, ncells * sizeof(int));

    // Loop over particles and record which grid cells contain particles
    for (int part = 0; part < npart; part++) {
        px = *(double *) ((char *) pos + part * s_pos[0] + 0 * s_pos[1]);
        py = *(double *) ((char *) pos + part * s_pos[0] + 1 * s_pos[1]);
        pz = *(double *) ((char *) pos + part * s_pos[0] + 2 * s_pos[1]);

        ix = floor(px / cellsize);
        iy = floor(py / cellsize);
        iz = floor(pz / cellsize);

        // Print a warning if a particle is given outside of the acceptable range:
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz) {
            std::cout << "part: " << part << " px, py, pz: " << px << "," << py << "," << pz << std::endl;
            std::cout << "ix, iy, iz: " << ix << "," << iy << "," << iz << std::endl;
            PyErr_SetString(PyExc_ValueError, "WARNING: Particle out of range [0, boxsize]");
        } else {
            grid_index = ((iz * ny + iy) * nx) + ix;
            grid[grid_index] += 1;
        }
    }

    int p = 0;
    auto return_grid = static_cast<double *>(malloc(ncells * sizeof(double) * 3));

    for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
            for (iz = 0; iz < nz; iz++) {
                grid_index = ((iz * ny + iy) * nx) + ix;

                // If cell unoccupied, add box particle
                if (grid[grid_index] == 0) {
                    return_grid[p * 3 + 0] = (ix + 0.5) * cellsize;
                    return_grid[p * 3 + 1] = (iy + 0.5) * cellsize;
                    return_grid[p * 3 + 2] = (iz + 0.5) * cellsize;

                    p++;
                }
            }

    // Convert the grid to a Numpy Array
    free(grid);
    return_grid = static_cast<double *>(realloc(return_grid, p * sizeof(double) * 3));
    auto res = create_numpy_array(return_grid, p, 3);
    free(return_grid);

    printf("Created %d particles to fill grid.\n", p);

    return res;

}

PyObject *convert_to_healpix_implementation(PyObject *dict_in, double boxsize, PyObject *centers_py,
                                            int randomizeshells, int randomizeradii, double pmass) {


    double num_base_pixels = 12;
    int npart_max = 1e8;

    //
    seed();

    auto r = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::RADIUS]);
    auto r_data = convert_to_std_vector(r);
    auto u = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::INTERNALENERGY]);
    auto u_data = convert_to_std_vector(u);
    auto pres = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::PRESSURE]);
    auto p_data = convert_to_std_vector(pres);
    auto rho = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::DENSITY]);
    auto rho_data = convert_to_std_vector(rho);
    auto mr = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::MR]);
    auto mr_data = convert_to_std_vector(mr);
    auto xnuc_radial = (PyArrayObject *) PyDict_GetItemString(dict_in, f[N::NUCLEARCOMPOSITION]);

    auto npoints = PyArray_DIMS(r)[0];
    int nspecies = PyArray_DIMS(xnuc_radial)[1];
//    auto centers = (double *) centers_py;
    double centers0 = PyFloat_AsDouble(PyList_GetItem(centers_py, 0));
    double centers1 = PyFloat_AsDouble(PyList_GetItem(centers_py, 1));
    double centers2 = PyFloat_AsDouble(PyList_GetItem(centers_py, 2));
    double centers[3] = {centers0, centers1, centers2};

    std::vector<double> xc_data, xo_data;
    auto s = PyArray_STRIDES(xnuc_radial);
    for (int i = 0; i < npoints; i++) {
        double xchere = * (double *) ((char *) PyArray_DATA(xnuc_radial) + i * s[0] + 1 * s[1]);
        double xohere = * (double *) ((char *) PyArray_DATA(xnuc_radial) + i * s[0] + 2 * s[1]);
        xc_data.push_back(xchere);
        xo_data.push_back(xohere);
    }

    tk::spline spline_carbon(r_data, xc_data);
    tk::spline spline_oxygen(r_data, xo_data);
    tk::spline spline_u(r_data, u_data);
    tk::spline spline_p(r_data, p_data);
    tk::spline spline_rho(r_data, rho_data);

    auto three_vector_data_size = 3 * npart_max * sizeof(double);
    auto scalar_data_size = npart_max * sizeof(double);
    auto ndata_pos = (double *) malloc(three_vector_data_size);
    auto ndata_vel = (double *) malloc(three_vector_data_size);
    auto ndata_pres = (double *) malloc(scalar_data_size);
    auto ndata_mass = (double *) malloc(scalar_data_size);
    auto ndata_rho = (double *) malloc(scalar_data_size);
    auto ndata_u = (double *) malloc(scalar_data_size);
    auto ndata_xnuc = (double *) malloc(nspecies * npart_max * sizeof(double));

    double shell_min_pmass = num_base_pixels * pmass;
    double rollover_mass = 0;
    double shell_radius;
    double shell_mass;
    double width = 0;
    double phis[3];
    double hp_vector[3];
    double i_prev = 0;
    int n_pix;
    long nside;
    double n2;

    int p = 0;


    for (int i = 0; i < (npoints - 1); i++) {

        if (i == 0) {
            shell_mass = mr_data[0];
        } else {
            shell_mass = mr_data[i + 1] - mr_data[i] + rollover_mass;
        }

        if (shell_mass > shell_min_pmass || i == npoints - 2) {

            if (i == npoints - 2) {
                shell_mass = shell_min_pmass;
            }

            double n_cells = shell_mass / pmass;
            n2 = sqrt(n_cells / num_base_pixels);
            nside = floor(n2 + 0.5);
            n_pix = static_cast<int>(num_base_pixels * pow(nside, 2));
            shell_radius = 0.5 * (r_data[i_prev] + r_data[i]);

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

                    double x = hp_vector[0];
                    double y = cos(phis[0]) * hp_vector[1] - sin(phis[0]) * hp_vector[2];
                    double z = sin(phis[0]) * hp_vector[1] + cos(phis[0]) * hp_vector[2];

                    double x2 = cos(phis[1]) * x + sin(phis[1]) * z;
                    double y2 = y;
                    double z2 = -sin(phis[1]) * x + cos(phis[1]) * z;

                    ndata_pos[p * 3 + 0] = shell_radius * (cos(phis[2]) * x2 - sin(phis[2]) * y2) + centers[0];
                    ndata_pos[p * 3 + 1] = shell_radius * (sin(phis[2]) * x2 + cos(phis[2]) * y2) + centers[1];
                    ndata_pos[p * 3 + 2] = shell_radius * z2 + centers[2];
                } else {
                    for (int k = 0; k < 3; k++)
                        ndata_pos[p * 3 + k] = hp_vector[k] * shell_radius + centers[k];
                }


                // Assign the values of the primitive variables to particle p
                ndata_mass[p] = shell_mass / n_pix;
                ndata_rho[p] = spline_rho(shell_radius);
                ndata_u[p] = spline_u(shell_radius);
                ndata_pres[p] = spline_p(shell_radius);
                double xCHere = spline_carbon(shell_radius);
                double xOHere = spline_oxygen(shell_radius);

                for (int k = 0; k < nspecies; k++) {
                    ndata_xnuc[p * nspecies + k] = 0;
                }

                if (xCHere + xOHere >= 1) { xOHere = 1 - xCHere; }
                ndata_xnuc[p * nspecies + 0] = 1 - xCHere - xOHere; // He
                ndata_xnuc[p * nspecies + 1] = xCHere; // C
                ndata_xnuc[p * nspecies + 2] = xOHere; // O

                // Set a density floor
                if (ndata_rho[p] < 1e-5) {
                    ndata_rho[p] = 1e-5;
                }

                // Iterate the particle counter
                p++;
            }

            rollover_mass = 0;
            i_prev = i;
        } else {
            rollover_mass = shell_mass;
        }
    }

    memset(ndata_vel, 0, 3 * p * sizeof(double));
    auto *particleIDs = static_cast<uint *>(malloc(p * sizeof(uint)));
    for (int l = 0; l < p; l++) { particleIDs[l] = l + 1; }

    printf("Created %d particles.\n", p);


    auto dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::COORDINATES], (PyObject *) create_numpy_array(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, f[N::MASSES], (PyObject *) create_numpy_array(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, f[N::DENSITY], (PyObject *) create_numpy_array(ndata_rho, p, 1));
    free(ndata_rho);
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject *) create_numpy_array(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject *) create_numpy_array(ndata_pres, p, 1));
    free(ndata_pres);
    PyDict_SetStolenItem(dict, f[N::VELOCITIES], (PyObject *) create_numpy_array(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, f[N::NUCLEARCOMPOSITION], (PyObject *) create_numpy_array(ndata_xnuc, p, nspecies));
    free(ndata_xnuc);

    return dict;
}

PyObject *convert_to_healpix(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyObject *dict;
    double boxsize;
    int randomizeshells = 0;
    int randomizeradii = 0;
    double pmass = 1e-6 * msol;
    PyObject *centers = nullptr;

    const char *kwlist[] = {"dict", "boxsize", "centers", "randomizeshells", "randomizeradii", "pmass", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "Od|Oiid:convert_to_healpix("
                                     "dict, boxsize, "
                                     "[centers, randomizeshells, randomizeradii, pmass])",
                                     keywords,
                                     &dict, &boxsize, &centers, &randomizeshells, &randomizeradii, &pmass)) {
        return nullptr;
    }

    auto ret_dict = convert_to_healpix_implementation(dict,
                                                      boxsize,
                                                      centers,
                                                      randomizeshells,
                                                      randomizeradii,
                                                      pmass);

    return ret_dict;
}

PyMODINIT_FUNC PyInit_create_ics(void) {
    import_array();

    return PyModule_Create(&moduledef_create_ics);
}
