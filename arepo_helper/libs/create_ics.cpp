#include <Python.h>
#include <arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "helm_eos.h"
#include "mersenne.h"
#include "pyhelm_eos.h"
#include "pix2vec_ring.h"
#include "utils.h"
#include "create_ics.h"
#include "spline.h"
#include "const.h"

int compare_points(const void *o1, const void *o2) {
    t_point *p1, *p2;
    p1 = (t_point *) o1;
    p2 = (t_point *) o2;

    if (p1->r < p2->r) return -1;
    else if (p1->r == p2->r) return 0;
    else return 1;
}

PyObject *create_particles_cube(PyObject *self, PyObject *args) {
    PyObject *data, *dict;
    t_helm_eos_table *helm_eos_table;
    long long bs, i, ix, iy, iz;
    int npart, usecells, ncells, nspecies;
    double temp, pmass, mass, dmass;
    PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
    double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
    t_point *cube;
    int j, k, n, index;
    double bsh;
    double *ndata_pos, *ndata_mass, *ndata_u, *ndata_vel, *ndata_xnuc;
    int p, cellspercell;
    double last_radius, next_radius, lr3, nr3, radius, newradius;
    struct eos_result res{};

    usecells = 0;
    nspecies = 3;
    temp = 0.0;
    pmass = 0.0;
    if (!PyArg_ParseTuple(args,
                          "OO&i|iidd:create_particles_cube( data, eos, npart, [usecells, nspecies, temp, pmass] )",
                          &data, &pyConvertHelmEos, &helm_eos_table, &npart, &usecells, &nspecies, &temp, &pmass)) {
        return nullptr;
    }

    data_u = data_H = data_HE = data_xnuc = nullptr;
    dm = (PyArrayObject *) PyDict_GetItemString(data, "dm");
    ncells = PyArray_DIMS(dm)[0];
    data_dm = (double *) PyArray_DATA(dm);
    r = (PyArrayObject *) PyDict_GetItemString(data, "r");
    data_r = (double *) PyArray_DATA(r);
    rho = (PyArrayObject *) PyDict_GetItemString(data, "rho");
    data_rho = (double *) PyArray_DATA(rho);

    if (temp == 0.0) {
        u = (PyArrayObject *) PyDict_GetItemString(data, "u");
        data_u = (double *) PyArray_DATA(u);
    }

    if (nspecies == 3) {
        H = (PyArrayObject *) PyDict_GetItemString(data, "H");
        data_H = (double *) PyArray_DATA(H);
        HE = (PyArrayObject *) PyDict_GetItemString(data, "HE");
        data_HE = (double *) PyArray_DATA(HE);
    } else {
        xnuc = (PyArrayObject *) PyDict_GetItemString(data, "xnuc");
        data_xnuc = (double *) PyArray_DATA(xnuc);
    }

    mass = 0.0;
    for (i = 0; i < ncells; i++) mass += data_dm[i];

    if (pmass > 0) {
        npart = floor(mass / pmass + 0.5);
    } else {
        pmass = mass / npart;
    }

    bs = floor(pow(3.0 * npart, 1. / 3.));
    bsh = (static_cast<double>(bs) - 1.0) / 2.0;
    cube = (t_point *) malloc(bs * bs * bs * sizeof(t_point));
    for (ix = 0; ix < bs; ix++) {
        for (iy = 0; iy < bs; iy++) {
            for (iz = 0; iz < bs; iz++) {
                i = (iz * bs + iy) * bs + ix;
                cube[i].x = static_cast<double>(ix) - bsh;
                cube[i].y = static_cast<double>(iy) - bsh;
                cube[i].z = static_cast<double>(iz) - bsh;
                cube[i].r = sqrt(cube[i].x * cube[i].x + cube[i].y * cube[i].y + cube[i].z * cube[i].z);
            }
        }
    }

    qsort(cube, bs * bs * bs, sizeof(t_point), compare_points);

    ndata_pos = (double *) malloc(3 * npart * 2 * sizeof(double));
    ndata_mass = (double *) malloc(npart * 2 * sizeof(double));
    ndata_u = (double *) malloc(npart * 2 * sizeof(double));
    ndata_vel = (double *) malloc(3 * npart * 2 * sizeof(double));
    ndata_xnuc = (double *) malloc(nspecies * npart * 2 * sizeof(double));

    p = 0;
    last_radius = 0.0;

    if (!usecells || usecells > ncells) {
        usecells = ncells;
    }

    cellspercell = std::max(1.0, floor(ncells / usecells + 0.5));
    if (usecells && cellspercell > 1) {
        usecells = floor(ncells / cellspercell);
        printf("Cells per cell: %d, using %d artificial cells.\n", cellspercell, usecells);
    }

    mass = 0;
    for (i = 0; i < usecells; i++) {
        for (j = static_cast<int>(i) * cellspercell; j < (i + 1) * cellspercell; j++)
            mass += data_dm[j];

        dmass = mass - p * pmass;
        n = floor(dmass / pmass + 0.5);

        next_radius = data_r[(i + 1) * cellspercell];
        lr3 = pow(last_radius, 3.0);
        nr3 = pow(next_radius, 3.0);

        index = static_cast<int>(i * cellspercell);
        for (j = 0; j < n; j++) {
            radius = cube[p].r;
            newradius = pow(1.0 * (j + 1.) / n * (nr3 - lr3) + lr3, 1.0 / 3.0);

            while (data_r[index] < newradius)
                index++;

            if (radius > 0) {
                ndata_pos[p * 3] = cube[p].x / radius * newradius;
                ndata_pos[p * 3 + 1] = cube[p].y / radius * newradius;
                ndata_pos[p * 3 + 2] = cube[p].z / radius * newradius;
            } else {
                ndata_pos[p * 3] = cube[p].x;
                ndata_pos[p * 3 + 1] = cube[p].y;
                ndata_pos[p * 3 + 2] = cube[p].z;
            }

            ndata_mass[p] = pmass;
            if (nspecies == 3) {
                ndata_xnuc[p * 3] = data_H[index];
                ndata_xnuc[p * 3 + 1] = data_HE[index];
                ndata_xnuc[p * 3 + 2] = 1. - data_H[index] - data_HE[index];
            } else {
                for (k = 0; k < nspecies; k++)
                    ndata_xnuc[p * nspecies + k] = data_xnuc[k];
            }

            if (temp > 0.0) {
                eos_calc_tgiven(helm_eos_table, data_rho[index], &ndata_xnuc[p], temp, &res);
                ndata_u[p] = res.e.v;
            } else {
                ndata_u[p] = data_u[index];
            }

            p++;
        }
        last_radius = next_radius;
    }

    memset(ndata_vel, 0, 3 * p * sizeof(double));

    free(cube);

    dict = PyDict_New();
    PyDict_SetStolenItem(dict, "ndata_pos", (PyObject *) createPyArray(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, "ndata_mass", (PyObject *) createPyArray(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, "ndata_u", (PyObject *) createPyArray(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, "ndata_vel", (PyObject *) createPyArray(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, "ndata_xnuc", (PyObject *) createPyArray(ndata_xnuc, p, nspecies));
    free(ndata_xnuc);
    PyDict_SetStolenItem(dict, "count", (PyObject *) PyLong_FromLong(p));

    return dict;
}

PyObject *create_particles_fill_grid(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyArrayObject *pyPos;
    double *pos, *grid_new;
    double boxsize, cellsize;
    double boxx, boxy, boxz;
    int *grid;
    int npart, part, pp, ncells, boxres;
    int ix, iy, iz, idx, nx, ny, nz;
    double px, py, pz;
    double longx, longy, longz;

    const char *kwlist[] = {"pos", "boxsize", "boxres", "longx", "longy", "longz", nullptr};
    auto keywords = (char **) kwlist;

    longx = 1.0;
    longy = 1.0;
    longz = 1.0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O!di|ddd:create_particles_fill_grid( pos, boxsize, boxres, [longx, longy, longz] )",
                                     keywords, &PyArray_Type, &pyPos, &boxsize, &boxres, &longx, &longy, &longz)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyPos) != 2 || PyArray_DIMS(pyPos)[1] != 3 || PyArray_TYPE(pyPos) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "pos has to be of dimensions [:,3] and type double");
        return nullptr;
    }

    pos = (double *) PyArray_DATA(pyPos);
    npart = PyArray_DIMS(pyPos)[0];

    printf("%d particles found.\n", npart);

    cellsize = boxsize / boxres;
    boxx = boxsize * longx;
    boxy = boxsize * longy;
    boxz = boxsize * longz;

    nx = floor(boxx / cellsize);
    ny = floor(boxy / cellsize);
    nz = floor(boxz / cellsize);

    printf("Building grid with %d x %d x %d cells.\n", nx, ny, nz);
    ncells = nx * ny * nz;

    grid = static_cast<int *>(malloc(ncells * sizeof(int)));
    memset(grid, 0, ncells * sizeof(int));

    for (part = 0; part < npart; part++) {
        px = *(double *) ((char *) pos + part * PyArray_STRIDES(pyPos)[0] + 0 * PyArray_STRIDES(pyPos)[1]);
        py = *(double *) ((char *) pos + part * PyArray_STRIDES(pyPos)[0] + 1 * PyArray_STRIDES(pyPos)[1]);
        pz = *(double *) ((char *) pos + part * PyArray_STRIDES(pyPos)[0] + 2 * PyArray_STRIDES(pyPos)[1]);

        ix = floor(px / cellsize);
        iy = floor(py / cellsize);
        iz = floor(pz / cellsize);

        if (ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz) {
            grid[((iz * ny + iy) * nx) + ix] += 1;
        }
    }

    pp = 0;
    grid_new = static_cast<double *>(malloc(ncells * sizeof(double) * 3));

    idx = 0;

    for (ix = 0; ix < nx; ix++)
        for (iy = 0; iy < ny; iy++)
            for (iz = 0; iz < nz; iz++) {
                idx = ((iz * ny + iy) * nx) + ix;
                if (grid[idx] == 0) {
                    grid_new[pp * 3] = (ix + 0.5) * cellsize;
                    grid_new[pp * 3 + 1] = (iy + 0.5) * cellsize;
                    grid_new[pp * 3 + 2] = (iz + 0.5) * cellsize;

                    pp++;
                }
            }

    free(grid);
    grid_new = static_cast<double *>(realloc(grid_new, pp * sizeof(double) * 3));

    printf("Created %d particles to fill grid.\n", pp);

    auto *res = (PyObject *) createPyArray(grid_new, pp, 3);
    free(grid_new);

    return res;
}

PyObject *convert_to_healpix_implementation(PyObject *wdec_dict,
                                            int nspecies,
                                            double boxsize,
                                            PyArrayObject *centers_py,
                                            int makebox,
                                            int randomizeshells,
                                            int randomizeradii,
                                            double pmass) {


    double num_base_pixels  = 12;
    double grid_energy      = 0;
    double grid_density     = 1e-4;
    double grid_temp        = 1500;
    int boxres              = 32;
    int npart_max           = 1e8;

    seed();

    auto r = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "r");
    auto r_data = convert_to_std_vector(r);
    auto u = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "u");
    auto u_data = convert_to_std_vector(u);
    auto rho = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "rho");
    auto rho_data = convert_to_std_vector(rho);
    auto temp = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "temp");
    auto temp_data = convert_to_std_vector(temp);
    auto mr = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "mr");
    auto mr_data = convert_to_std_vector(mr);
    auto xc = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "xc");
    auto xc_data = convert_to_std_vector(xc);
    auto xo = (PyArrayObject *) PyDict_GetItemString(wdec_dict, "xo");
    auto xo_data = convert_to_std_vector(xo);

    auto npoints = PyArray_DIMS(xc)[0];
    auto centers = (double *) PyArray_DATA(centers_py);

    tk::spline spline_carbon(r_data, xc_data);
    tk::spline spline_oxygen(r_data, xo_data);
    tk::spline spline_u(r_data, u_data);
    tk::spline spline_rho(r_data, rho_data);
    tk::spline spline_temp(r_data, temp_data);

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
    int min_nside = 1;
    double inner_radius, outer_radius, shell_radius;
    double interior_mass;
    double width = 0;
    double phis[3], hp_vector[3];
    int n_pix;
    long nside;
    double n2;

    int p = 0;

    for (int i = 0; i < (npoints - 1); i++) {


        inner_radius = r_data[i];
        outer_radius = r_data[i + 1];
        shell_radius = 0.5 * (inner_radius + outer_radius);
        interior_mass = mr_data[i];

        // TODO pmass might not be the best way to achieve this
        n2      = sqrt(interior_mass / pmass / num_base_pixels); /* increases with index */
        nside   = floor(n2 + 0.5);
        if (nside >= max_nside) { nside = max_nside;}
        if (nside < min_nside) {nside = min_nside;}
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

    auto dict = PyDict_New();
    PyDict_SetStolenItem(dict, "ndata_pos", (PyObject *) createPyArray(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, "ndata_mass", (PyObject *) createPyArray(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, "ndata_rho", (PyObject *) createPyArray(ndata_rho, p, 1));
    free(ndata_rho);
    PyDict_SetStolenItem(dict, "ndata_u", (PyObject *) createPyArray(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, "ndata_vel", (PyObject *) createPyArray(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, "ndata_temp", (PyObject *) createPyArray(ndata_temp, p, 1));
    free(ndata_temp);
    PyDict_SetStolenItem(dict, "ndata_xnuc", (PyObject *) createPyArray(ndata_xnuc, p, nspecies));
    free(ndata_xnuc);
    PyDict_SetStolenItem(dict, "count", (PyObject *) PyLong_FromLong(p));
    PyDict_SetStolenItem(dict, "boxsize", PyFloat_FromDouble(boxsize));

    return dict;
}

PyObject *convert_to_healpix(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyObject *wdec_dict;
    int nspecies;
    double boxsize;
    int makebox         = 1;
    int randomizeshells = 0;
    int randomizeradii  = 0;
    double pmass        = 1e-6 * msol;
    PyArrayObject *centers = nullptr;

    const char *kwlist[] = {"wdec_dict", "nspecies", "boxsize",
                            "centers", "makebox", "randomizeshells", "randomizeradii", "pmass", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "Oid|O!iiid:convert_to_healpix( "
                                     "wdec_dict, nspecies, boxsize, "
                                     "[centers, makebox, randomizeshells, randomizeradii, pmass] )",
                                     keywords,
                                     &wdec_dict, &nspecies, &boxsize,
                                     &PyArray_Type, &centers,
                                     &makebox, &randomizeshells, &randomizeradii, &pmass)) {
        return nullptr;
    }

    if (centers == nullptr) {
        auto cs = (double *) PyArray_DATA(centers);
        cs[0] = boxsize / 2;
        cs[1] = boxsize / 2;
        cs[2] = boxsize / 2;
    }

    auto dict = convert_to_healpix_implementation(wdec_dict,
                                                  nspecies,
                                                  boxsize,
                                                  centers,
                                                  makebox,
                                                  randomizeshells,
                                                  randomizeradii,
                                                  pmass);

    return dict;
}

// Python Module Definition
static PyMethodDef create_ics_Methods[] = {
        {"create_particles_cube", create_particles_cube,
        METH_VARARGS,""
         },
        {"convert_to_healpix", (PyCFunction) convert_to_healpix,
        METH_VARARGS | METH_KEYWORDS, ""
        },
        {"create_particles_fill_grid", (PyCFunction) create_particles_fill_grid,
         METH_VARARGS | METH_KEYWORDS,""
         },
        {nullptr, nullptr,
         0,nullptr
        }
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "create_ics",
        nullptr,
        -1,
        create_ics_Methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
};

PyMODINIT_FUNC PyInit_create_ics(void) {
    import_array();

    return PyModule_Create(&moduledef);
}
