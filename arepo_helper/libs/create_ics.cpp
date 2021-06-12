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
    PyDict_SetStolenItem(dict, "pos", (PyObject *) createPyArray(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, "mass", (PyObject *) createPyArray(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, "u", (PyObject *) createPyArray(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, "vel", (PyObject *) createPyArray(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, "xnuc", (PyObject *) createPyArray(ndata_xnuc, p, nspecies));
    free(ndata_xnuc);
    PyDict_SetStolenItem(dict, "count", (PyObject *) PyLong_FromLong(p));

    return dict;
}

PyObject *create_particles_healpix(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *data, *dict;
    t_helm_eos_table *helm_eos_table;
    PyArrayObject *dm, *r, *u, *xnuc, *H, *HE, *rho;
    double *data_dm, *data_r, *data_u, *data_xnuc, *data_H, *data_HE, *data_rho;
    double pmass, mass, mtot, masslast, width;
    int npart, nspecies, ncells, np;
    double cx, cy, cz;
    double *ndata_pos, *ndata_mass, *ndata_rho, *ndata_u, *ndata_vel, *ndata_xnuc;
    int idxlast, i, j, k, p, found, index, npart_allocated;
    double n1, n2, n, nlast, rad1, rad2, rad, vec[3];
    int makebox, boxres, nbox, ix, iy, iz;
    double boxsize, hboxsize, boxcellsize, minenergy;
    int randomizeshells, randomizeradii, noxnuc, fixgridpressure, transition_done;
    double phi1, phi2, phi3, x, y, z, x2, y2, z2;
    double gridenergy, griddensity, gridtemp, radius, boxfactor;
    double radius_last, dr, vol, volboxcell, volcell;
    double transition_radius, transition_pmass;

    const char *kwlist[] = {"data", "eos", "npart", "nspecies", "pmass", "makebox", "boxsize", "boxres",
                             "randomizeshells", "randomizeradii", "minenergy", "fixgridpressure", "griddensity",
                             "boxfactor", "gridtemp", "cx", "cy", "cz", "transition_radius", "transition_pmass", nullptr};
    auto keywords = (char **) kwlist;

    nspecies = 3;
    pmass = 0;
    makebox = 0;
    boxsize = 0;
    boxres = 16;
    randomizeshells = 0;
    randomizeradii = 0;
    minenergy = 0;
    fixgridpressure = 0;
    griddensity = 1e-5;
    boxfactor = 10;
    gridtemp = 1e8;
    cx = cy = cz = 0;
    transition_radius = 0;
    transition_pmass = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "OO&i|ididiiididddddddd:create_particles_healpix( data, eos, npart, [nspecies, pmass, makebox, boxsize, boxres, randomizeshells, randomizeradii, minenergy, fixgridpressure, griddensity, boxfactor, gridtemp, cx, cy, cz, transition_radius, transition_pmass] )",
                                     keywords, &data, &pyConvertHelmEos, &helm_eos_table, &npart, &nspecies, &pmass,
                                     &makebox, &boxsize, &boxres, &randomizeshells, &randomizeradii, &minenergy,
                                     &fixgridpressure, &griddensity, &boxfactor, &gridtemp, &cx, &cy, &cz,
                                     &transition_radius, &transition_pmass)) {
        return nullptr;
    }

    dm = (PyArrayObject *) PyDict_GetItemString(data, "dm");
    ncells = PyArray_DIMS(dm)[0];
    data_dm = (double *) PyArray_DATA(dm);
    r = (PyArrayObject *) PyDict_GetItemString(data, "r");
    data_r = (double *) PyArray_DATA(r);
    rho = (PyArrayObject *) PyDict_GetItemString(data, "rho");
    data_rho = (double *) PyArray_DATA(rho);
    u = (PyArrayObject *) PyDict_GetItemString(data, "u");
    data_u = (double *) PyArray_DATA(u);

    data_H = data_HE = data_xnuc = nullptr;
    if (nspecies == 3 && !PyDict_Contains(data, PyUnicode_FromString("xnuc"))) {
        noxnuc = 1;
        H = (PyArrayObject *) PyDict_GetItemString(data, "H");
        data_H = (double *) PyArray_DATA(H);
        HE = (PyArrayObject *) PyDict_GetItemString(data, "HE");
        data_HE = (double *) PyArray_DATA(HE);
    } else {
        noxnuc = 0;
        xnuc = (PyArrayObject *) PyDict_GetItemString(data, "xnuc");
        data_xnuc = (double *) PyArray_DATA(xnuc);
    }

    mtot = 0.0;
    for (i = 0; i < ncells; i++) mtot += data_dm[i];

    printf("Total mass: %g solar masses\n", mtot / 1.989e33);

    if (pmass > 0) {
        npart = floor(mtot / pmass + 0.5);
    } else {
        pmass = mtot / npart;
    }

    if (transition_radius > 0) {
        transition_done = 0;
        printf("Doing transition at r=%g to pmass=%g\n", transition_radius, transition_pmass);
    } else {
        transition_done = 1;
        printf("Using %d particles of mass %g (%g solar masses)\n", npart, pmass, pmass / 1.989e33);
    }

    seed();

    npart_allocated = npart;
    ndata_pos = (double *) malloc(3 * npart_allocated * sizeof(double));
    ndata_mass = (double *) malloc(npart_allocated * sizeof(double));
    ndata_rho = (double *) malloc(npart_allocated * sizeof(double));
    ndata_u = (double *) malloc(npart_allocated * sizeof(double));
    ndata_vel = (double *) malloc(3 * npart_allocated * sizeof(double));
    ndata_xnuc = (double *) malloc(nspecies * npart_allocated * sizeof(double));
    p = 0;

    idxlast = 0;
    nlast = 0;
    masslast = 0;
    mass = 0;
    rad1 = data_r[0];
    rad2 = rad1;
    found = 0;
    volcell = 0;
    phi1 = phi2 = phi3 = 0;
    np = 0;

    i = 0;
    while (i < ncells - 1) {
        rad2 = data_r[i + 1]; /* increases with index */
        n1 = sqrt(M_PI / 12.) * (rad2 + rad1) / (rad2 - rad1); /* decreases with index */

        if ((!transition_done) && (rad2 > transition_radius)) {
            pmass = transition_pmass;
            transition_done = 1;
            printf("Transition done at rad2=%g, index=%d, mass=%g, particles before: %d.\n", rad2, i, mass / 1.989e33,
                   p);
        }

        mass += data_dm[i]; /* increases with index */
        n2 = sqrt(mass / pmass / 12.); /* increases with index */

        if (floor(n2) > nlast) {
            nlast = floor(n2);
            idxlast = i;
            masslast = mass;
        }

        if (n2 > n1 && !found) {
            n = floor(n2 + 0.5);
            found = 1;
        }

        if (found && nlast >= n) {
            i = idxlast;
            mass = masslast;
            rad2 = data_r[i + 1];

            rad = 0.5 * (rad2 + rad1);
            width = rad2 - rad1;
            index = i;
            while (data_r[index] > rad) index--;

            np = static_cast<int>(12 * n * n);

            if (randomizeshells) {
                phi1 = randMT() * 2. * M_PI;
                phi2 = randMT() * 2. * M_PI;
                phi3 = randMT() * 2. * M_PI;
            }

            while (p + np >= npart_allocated) {
                npart_allocated *= 2;

                ndata_pos = (double *) realloc(ndata_pos, 3 * npart_allocated * sizeof(double));
                ndata_mass = (double *) realloc(ndata_mass, npart_allocated * sizeof(double));
                ndata_rho = (double *) realloc(ndata_rho, npart_allocated * sizeof(double));
                ndata_u = (double *) realloc(ndata_u, npart_allocated * sizeof(double));
                ndata_vel = (double *) realloc(ndata_vel, 3 * npart_allocated * sizeof(double));
                ndata_xnuc = (double *) realloc(ndata_xnuc, nspecies * npart_allocated * sizeof(double));
            }

            for (j = 0; j < np; j++) {
                auto n_int = (int) n;
                pix2vec_ring(n_int, j, vec);

                double prad = rad;
                if (randomizeradii)
                    prad += 0.1 * width * (randMT() - 0.5);

                if (randomizeshells) {
                    x = vec[0];
                    y = cos(phi1) * vec[1] - sin(phi1) * vec[2];
                    z = sin(phi1) * vec[1] + cos(phi1) * vec[2];

                    x2 = cos(phi2) * x + sin(phi2) * z;
                    y2 = y;
                    z2 = -sin(phi2) * x + cos(phi2) * z;

                    ndata_pos[p * 3 + 0] = prad * (cos(phi3) * x2 - sin(phi3) * y2);
                    ndata_pos[p * 3 + 1] = prad * (sin(phi3) * x2 + cos(phi3) * y2);
                    ndata_pos[p * 3 + 2] = prad * z2;
                } else {
                    for (k = 0; k < 3; k++)
                        ndata_pos[p * 3 + k] = vec[k] * prad;
                }

                if (makebox) {
                    int inbox = 1;
                    for (k = 0; k < 3; k++)
                        if (ndata_pos[p * 3 + k] < -0.5 * boxsize || ndata_pos[p * 3 + k] > +0.5 * boxsize) {
                            inbox = 0;
                            break;
                        }
                    if (!inbox)
                        continue;
                }

                ndata_mass[p] = pmass;
                ndata_rho[p] = data_rho[index];
                ndata_u[p] = data_u[index];
                if (noxnuc) {
                    ndata_xnuc[p * 3] = data_H[index];
                    ndata_xnuc[p * 3 + 1] = data_HE[index];
                    ndata_xnuc[p * 3 + 2] = 1. - data_H[index] - data_HE[index];
                } else {
                    for (k = 0; k < nspecies; k++)
                        ndata_xnuc[p * nspecies + k] = data_xnuc[k];
                }

                p++;
            }

            volcell = 4. / 3. * M_PI * (rad2 * rad2 * rad2 - rad1 * rad1 * rad1) / np;

            rad1 = rad2;
            mass -= pmass * 12. * n * n;
            found = 0;
            nlast = 0;
        }

        i++;
    }

    n = floor(sqrt(mass / pmass / 12.) + 0.5);
    /* put last point near radius of star
    if (ncells > 50) {
      rad1 = data_r[i - 6];
    } else {
      rad1 = data_r[i - 1];
    }
    */
    rad = 0.5 * (rad2 + rad1);
    width = rad2 - rad1;
    index = i;
    while (data_r[index] > rad) index--;

    if (randomizeshells) {
        phi1 = randMT() * 2. * M_PI;
        phi2 = randMT() * 2. * M_PI;
        phi3 = randMT() * 2. * M_PI;
    }

    double prad = rad;
    if (randomizeradii)
        prad += 0.1 * width * (randMT() - 0.5);

    while (p + np >= npart_allocated) {
        npart_allocated *= 2;

        ndata_pos = (double *) realloc(ndata_pos, 3 * npart_allocated * sizeof(double));
        ndata_mass = (double *) realloc(ndata_mass, npart_allocated * sizeof(double));
        ndata_rho = (double *) realloc(ndata_rho, npart_allocated * sizeof(double));
        ndata_u = (double *) realloc(ndata_u, npart_allocated * sizeof(double));
        ndata_vel = (double *) realloc(ndata_vel, 3 * npart_allocated * sizeof(double));
        ndata_xnuc = (double *) realloc(ndata_xnuc, nspecies * npart_allocated * sizeof(double));
    }

    np = static_cast<int>(12 * n * n);
    for (j = 0; j < np; j++) {
        auto n_int = (int) n;
        pix2vec_ring(n_int, j, vec);

        if (randomizeshells) {
            x = vec[0];
            y = cos(phi1) * vec[1] - sin(phi1) * vec[2];
            z = sin(phi1) * vec[1] + cos(phi1) * vec[2];

            x2 = cos(phi2) * x + sin(phi2) * z;
            y2 = y;
            z2 = -sin(phi2) * x + cos(phi2) * z;

            ndata_pos[p * 3 + 0] = prad * (cos(phi3) * x2 - sin(phi3) * y2);
            ndata_pos[p * 3 + 1] = prad * (sin(phi3) * x2 + cos(phi3) * y2);
            ndata_pos[p * 3 + 2] = prad * z2;
        } else {
            for (k = 0; k < 3; k++)
                ndata_pos[p * 3 + k] = vec[k] * prad;
        }

        if (makebox) {
            int inbox = 1;
            for (k = 0; k < 3; k++)
                if (ndata_pos[p * 3 + k] < -0.5 * boxsize || ndata_pos[p * 3 + k] > +0.5 * boxsize) {
                    inbox = 0;
                    break;
                }
            if (inbox == 0)
                continue;
        }

        ndata_mass[p] = pmass;
        ndata_rho[p] = data_rho[index];

        ndata_u[p] = data_u[index];
        if (ndata_u[p] < minenergy) ndata_u[p] = minenergy;
        if (noxnuc) {
            ndata_xnuc[p * 3] = data_H[index];
            ndata_xnuc[p * 3 + 1] = data_HE[index];
            ndata_xnuc[p * 3 + 2] = 1. - data_H[index] - data_HE[index];
        } else {
            for (k = 0; k < nspecies; k++)
                ndata_xnuc[p * nspecies + k] = data_xnuc[k];
        }

        p++;
    }

    if (makebox) {
        auto *grid_xnuc = static_cast<double *>(malloc(nspecies * sizeof(double)));
        if (noxnuc) {
            grid_xnuc[0] = data_H[index];
            grid_xnuc[1] = data_HE[index];
            grid_xnuc[2] = 1. - data_H[index] - data_HE[index];
        } else {
            for (k = 0; k < nspecies; k++)
                grid_xnuc[k] = data_xnuc[k];
        }

        if (fixgridpressure) {
            struct eos_result res {};
            double pressure, temp;

            temp = -1.;
            eos_calc_egiven(helm_eos_table, data_rho[index], grid_xnuc, data_u[index], &temp, &res);
            printf("rho=%g, u=%g, temp=%g, p=%g\n", data_rho[index], data_u[index], temp, res.p.v);

            temp = -1;
            pressure = res.p.v;
            eos_calc_pgiven(helm_eos_table, griddensity, grid_xnuc, pressure, &temp, &res);
            printf("rho=%g, u=%g, temp=%g, p=%g\n", griddensity, res.e.v, temp, res.p.v);
            gridenergy = res.e.v;

            if (temp == -1)
                temp = gridtemp;
            eos_calc_tgiven(helm_eos_table, griddensity, grid_xnuc, temp, &res);
            printf("rho=%g, u=%g, temp=%g, p=%g\n", griddensity, res.e.v, temp, res.p.v);
            gridenergy = res.e.v;
        } else {
            gridenergy = minenergy;
        }

        hboxsize = boxsize / 2.;
        boxcellsize = boxsize / (double) boxres;
        volboxcell = boxcellsize * boxcellsize * boxcellsize;

        radius_last = rad2;
        //volcell = 4./3. * M_PI * (rad2*rad2*rad2 - rad1*rad1*rad1);
        printf("volcell=%g, volboxcell=%g\n", volcell, volboxcell);
        volcell *= boxfactor;
        while (volcell < volboxcell && boxfactor > 0.) {
            dr = pow(volcell, 1. / 3.);

            vol = 4. / 3. * M_PI * (pow(radius_last + dr, 3) - pow(radius_last, 3));
            n = floor(sqrt(vol / (12. * volcell)) + 0.5);
            np = static_cast<int>(12 * n * n);
            radius = radius_last + 0.5 * dr;

            printf("volcell=%g, volboxcell=%g, np=%d\n", volcell, volboxcell, np);

            while (p + np >= npart_allocated) {
                npart_allocated *= 2;

                ndata_pos = (double *) realloc(ndata_pos, 3 * npart_allocated * sizeof(double));
                ndata_mass = (double *) realloc(ndata_mass, npart_allocated * sizeof(double));
                ndata_rho = (double *) realloc(ndata_rho, npart_allocated * sizeof(double));
                ndata_u = (double *) realloc(ndata_u, npart_allocated * sizeof(double));
                ndata_vel = (double *) realloc(ndata_vel, 3 * npart_allocated * sizeof(double));
                ndata_xnuc = (double *) realloc(ndata_xnuc, nspecies * npart_allocated * sizeof(double));
            }

            if (randomizeshells) {
                phi1 = randMT() * 2. * M_PI;
                phi2 = randMT() * 2. * M_PI;
                phi3 = randMT() * 2. * M_PI;
            }

            for (j = 0; j < np; j++) {
                auto n_int = (int) n;
                pix2vec_ring(n_int, j, vec);

                if (randomizeshells) {
                    x = vec[0];
                    y = cos(phi1) * vec[1] - sin(phi1) * vec[2];
                    z = sin(phi1) * vec[1] + cos(phi1) * vec[2];

                    x2 = cos(phi2) * x + sin(phi2) * z;
                    y2 = y;
                    z2 = -sin(phi2) * x + cos(phi2) * z;

                    ndata_pos[p * 3 + 0] = radius * (cos(phi3) * x2 - sin(phi3) * y2);
                    ndata_pos[p * 3 + 1] = radius * (sin(phi3) * x2 + cos(phi3) * y2);
                    ndata_pos[p * 3 + 2] = radius * z2;
                } else {
                    for (k = 0; k < 3; k++)
                        ndata_pos[p * 3 + k] = vec[k] * radius;
                }

                int inbox = 1;
                for (k = 0; k < 3; k++)
                    if (ndata_pos[p * 3 + k] < -hboxsize || ndata_pos[p * 3 + k] > +hboxsize) {
                        inbox = 0;
                        break;
                    }
                if (inbox == 0)
                    continue;

                ndata_rho[p] = griddensity;
                ndata_mass[p] = griddensity * volcell;

                for (k = 0; k < nspecies; k++)
                    ndata_xnuc[p * nspecies + k] = grid_xnuc[k];

                p++;
            }

            radius_last += dr;
            volcell *= boxfactor;
        }

        for (i = 0; i < p; i++) {
            ndata_pos[i * 3 + 0] += hboxsize + cx;
            ndata_pos[i * 3 + 1] += hboxsize + cy;
            ndata_pos[i * 3 + 2] += hboxsize + cz;
        }

        nbox = boxres * boxres * boxres;

        while (p + nbox >= npart_allocated) {
            npart_allocated *= 2;

            ndata_pos = (double *) realloc(ndata_pos, 3 * npart_allocated * sizeof(double));
            ndata_mass = (double *) realloc(ndata_mass, npart_allocated * sizeof(double));
            ndata_rho = (double *) realloc(ndata_rho, npart_allocated * sizeof(double));
            ndata_u = (double *) realloc(ndata_u, npart_allocated * sizeof(double));
            ndata_vel = (double *) realloc(ndata_vel, 3 * npart_allocated * sizeof(double));
            ndata_xnuc = (double *) realloc(ndata_xnuc, nspecies * npart_allocated * sizeof(double));
        }

        char *box = static_cast<char *>(malloc(nbox));
        memset(box, 0, nbox);
        for (i = 0; i < p; i++) {
            ix = floor(ndata_pos[i * 3 + 0] / boxcellsize);
            iy = floor(ndata_pos[i * 3 + 1] / boxcellsize);
            iz = floor(ndata_pos[i * 3 + 2] / boxcellsize);

            box[iz * boxres * boxres + iy * boxres + ix] = 1;
        }

        int boxCount = 0;
        for (i = 0; i < nbox; i++) {
            ix = i % boxres;
            iy = (i / boxres) % boxres;
            iz = i / (boxres * boxres);

            x = (ix + 0.5) * boxcellsize;
            y = (iy + 0.5) * boxcellsize;
            z = (iz + 0.5) * boxcellsize;

            if (!box[i]) {
                ndata_pos[p * 3 + 0] = x;
                ndata_pos[p * 3 + 1] = y;
                ndata_pos[p * 3 + 2] = z;

                ndata_rho[p] = griddensity;
                ndata_mass[p] = griddensity * volboxcell;

                ndata_u[p] = gridenergy;
                for (k = 0; k < nspecies; k++)
                    ndata_xnuc[p * nspecies + k] = grid_xnuc[k];

                p++;
                boxCount++;
            }
        }
        free(box);
        printf("Added %d box cells.\n", boxCount);

        free(grid_xnuc);
    }

    memset(ndata_vel, 0, 3 * p * sizeof(double));

    printf("Created %d particles.\n", p);

    dict = PyDict_New();
    PyDict_SetStolenItem(dict, "pos", (PyObject *) createPyArray(ndata_pos, p, 3));
    free(ndata_pos);
    PyDict_SetStolenItem(dict, "mass", (PyObject *) createPyArray(ndata_mass, p, 1));
    free(ndata_mass);
    PyDict_SetStolenItem(dict, "rho", (PyObject *) createPyArray(ndata_rho, p, 1));
    free(ndata_rho);
    PyDict_SetStolenItem(dict, "u", (PyObject *) createPyArray(ndata_u, p, 1));
    free(ndata_u);
    PyDict_SetStolenItem(dict, "vel", (PyObject *) createPyArray(ndata_vel, p, 3));
    free(ndata_vel);
    PyDict_SetStolenItem(dict, "xnuc", (PyObject *) createPyArray(ndata_xnuc, p, nspecies));
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


// Python Module Definition
static PyMethodDef create_ics_Methods[] = {
        {"create_particles_cube",
                create_particles_cube,
         METH_VARARGS,
         ""
         },
        {"create_particles_healpix",
         (PyCFunction) create_particles_healpix,
         METH_VARARGS | METH_KEYWORDS,
         ""
         },
        {"create_particles_fill_grid",
         (PyCFunction) create_particles_fill_grid,
         METH_VARARGS | METH_KEYWORDS,
         ""
         },
        {nullptr,
         nullptr,
         0,
         nullptr
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
