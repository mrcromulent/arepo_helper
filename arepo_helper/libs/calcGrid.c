#include <Python.h>
#include <arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "omp_util.h"
#include "sph.h"

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define sqr(x) ((x)*(x))

inline void PyDict_SetStolenItem(PyObject *dict, const char *key, PyObject *object) {
        PyDict_SetItemString(dict, key, object);
        Py_DECREF(object);
}

inline double _getkernel( double h, double r2 ) {
        double coeff1, coeff2, coeff5;
        double hinv, hinv3, u;
        coeff1 = 8.0 / M_PI;
        coeff2 = coeff1 * 6.0;
        coeff5 = coeff1 * 2.0;

        hinv = 1.0 / h;
        hinv3 = hinv*hinv*hinv;
        u = sqrt(r2)*hinv;
        if (u < 0.5) {
                return hinv3 * ( coeff1 + coeff2*(u-1.0)*u*u );
        } else {
                return hinv3 * coeff5 * (1.0-u) * (1.0-u) * (1.0-u);
        }
}

inline double _getkernelewald( double h2, double r2 ) {
        double tmp = 1. - r2 / h2;
        return tmp * tmp;
}

PyObject* _calcGrid(PyObject *self, PyObject *args, PyObject *kwargs) {
        PyArrayObject *pos, *hsml, *mass, *rho, *value, *pyGrid;
        int npart, nx, ny, nz, cells, numthreads;
        int dims[3];
        double bx, by, bz, cx, cy, cz;
        double *grid;
        int part, proj, norm, ewaldkernel;
        double cellsizex, cellsizey, cellsizez;
        double start;
        char *kwlist[] = {"pos", "hsml", "mass", "rho", "value", "nx", "ny", "nz", "boxx", "boxy", "boxz", "centerx", "centery", "centerz", "proj", "norm", "ewaldkernel", "numthreads", NULL};
        
        start = get_time();

        proj = 0;
        norm = 0;
        ewaldkernel = 0;
        numthreads = 1;
        if (!PyArg_ParseTupleAndKeywords( args, kwargs, "O!O!O!O!O!iiidddddd|iiii:calcGrid( pos, hsml, mass, rho, value, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz, [proj, norm, ewaldkernel, numthreads] )", kwlist, &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &rho, &PyArray_Type, &value, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz, &proj, &norm, &ewaldkernel, &numthreads )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(rho) != 1 || PyArray_TYPE(rho) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "rho has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(rho)[0] || npart != PyArray_DIMS(value)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml, rho and value have to have the same size in the first dimension" );
                return 0;
        }


        if (npart == 0) {
          PyErr_SetString( PyExc_ValueError, "Number of particles has to be larger than zero" );
          return 0;
        }


        if (proj) {
                dims[0] = nx;
                dims[1] = ny;
                pyGrid = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
                grid = (double*)PyArray_DATA(pyGrid);
                cells = nx*ny;
        } else {
                dims[0] = nx;
                dims[1] = ny;
                dims[2] = nz;
                pyGrid = (PyArrayObject *)PyArray_FromDims( 3, dims, NPY_DOUBLE );
                grid = (double*)PyArray_DATA(pyGrid);
                cells = nx*ny*nz;               
        }
        memset( grid, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;
        cellsizez = bz / nz;

        set_num_threads(numthreads);

        printf("Mapping onto grid with %d thread(s)\n", numthreads);

        #pragma omp parallel private(part)
        {
                double px, py, pz, h, h2, m, r, v, cpx, cpy, cpz, r2, sum;
                int x, y, z0, z1;
                int xmin, xmax, ymin, ymax, zmin, zmax, zmid;
                double *data_pos, *data_hsml, *data_mass, *data_rho, *data_value;

                #pragma omp for schedule(dynamic, 1) nowait
                for (part=0; part<npart; part++) {
                        data_pos = (double*)((char*)PyArray_DATA(pos) + part*PyArray_STRIDES(pos)[0]);
                        px = *data_pos;
                        data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                        py = *data_pos;
                        data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                        pz = *data_pos;
                
                        data_hsml = (double*)((char*)PyArray_DATA(hsml) + part*PyArray_STRIDES(hsml)[0]);
                        h = *data_hsml;
                        h2 = h*h;

                        data_mass = (double*)((char*)PyArray_DATA(mass) + part*PyArray_STRIDES(mass)[0]);
                        m = *data_mass;

                        data_rho = (double*)((char*)PyArray_DATA(rho) + part*PyArray_STRIDES(rho)[0]);
                        r = *data_rho;

                        data_value = (double*)((char*)PyArray_DATA(value) + part*PyArray_STRIDES(value)[0]);
                        v = *data_value;

                        xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                        xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                        ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                        ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                        zmin = max( floor( (pz - h - cz + 0.5*bz) / cellsizez ), 0 );
                        zmax = min( ceil( (pz + h - cz + 0.5*bz) / cellsizez ), nz-1 );

                        zmid = floor( 0.5 * (zmin+zmax) + 0.5 );
        
                        if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && zmin < nz && zmax >= 0) {
                                if (norm) {
                                        sum = 0.;
                                        for (x=xmin; x<=xmax; x++) {
                                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                                for (y=ymin; y<=ymax; y++) {
                                                        cpy = -0.5*by + by*(y+0.5)/ny;
        
                                                        if (proj) {
                                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) );
                                                                if (r2 < h2) {
                                                                        if (ewaldkernel) {
                                                                                sum += h * _getkernelewald( h2, r2 );
                                                                        } else {
                                                                                sum += h * _getkernel( h, r2 );
                                                                        }
                                                                }
                                                        } else {                                
                                                                for (z0=zmid; z0>=zmin; z0--) {
                                                                        cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                                        r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                                        if (r2 > h2) break;
                                                                        sum += _getkernel( h, r2 );
                                                                }

                                                                for (z1=zmid+1; z1<=zmax; z1++) {
                                                                        cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                                        r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                                        if (r2 > h2) break;
                                                                        sum += _getkernel( h, r2 );
                                                                }
                                                        }
                                                }
                                        }
                                } else {
                                        sum = 1.0;
                                }
                        
                                for (x=xmin; x<=xmax; x++) {
                                        cpx = -0.5*bx + bx*(x+0.5)/nx;
                                        for (y=ymin; y<=ymax; y++) {
                                                cpy = -0.5*by + by*(y+0.5)/ny;
                                        
                                                if (proj) {
                                                        r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) );
                                                        if (r2 < h2) {
                                                                if (ewaldkernel) {
                                                                        #pragma omp atomic
                                                                        grid[x*ny + y] += h * _getkernelewald( h2, r2 ) * m * v / r / sum;
                                                                } else {
                                                                        #pragma omp atomic
                                                                        grid[x*ny + y] += h * _getkernel( h, r2 ) * m * v / r / sum;
                                                                }
                                                        }
                                                } else {                                
                                                        for (z0=zmid; z0>=zmin; z0--) {
                                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                                if (r2 > h2) break;
                                                                #pragma omp atomic
                                                                grid[(x*ny + y)*nz + z0] += _getkernel( h, r2 ) * m * v / r / sum;
                                                        }

                                                        for (z1=zmid+1; z1<=zmax; z1++) {
                                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                                if (r2 > h2) break;
                                                                #pragma omp atomic
                                                                grid[(x*ny + y)*nz + z1] += _getkernel( h, r2 ) * m * v / r / sum;
                                                        }
                                                }
                                        }
                                }       
                        }
                }
        }

        printf( "Calculation took %gs\n", get_time() - start );
        return PyArray_Return( pyGrid );
}

PyObject* _calcAbundGrid(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *abund, *pyGrid;
        int npart, nx, ny, nz, cells, cell;
        int dims[4];
        double bx, by, bz, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass, *data_abund, *xnuc;
        double *grid;
        int part, nspecies, sum, lost;
        double px, py, pz, h, h2, m, cpx, cpy, cpz, r2, mlost, mtot;
        int x, y, z0, z1, i;
        int xmin, xmax, ymin, ymax, zmin, zmax, zmid;
        double cellsizex, cellsizey, cellsizez, mincellsize;
        clock_t start;
        
        start = clock();

        if (!PyArg_ParseTuple( args, "O!O!O!O!iiidddddd:calcGrid( pos, hsml, mass, abund, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &abund, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [npart,3] and type double" );
                return 0;
        }
        npart = PyArray_DIMS(pos)[0];

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [npart] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [npart] and type double" );
                return 0;
        }

        if (PyArray_NDIM(abund) != 2 || PyArray_TYPE(abund) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "rho has to be of dimension [npart,nspecies] and type double" );
                return 0;
        }
        nspecies = PyArray_DIMS(abund)[1];

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(abund)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsmlm, mass and abund have to have the same size in the first dimension" );
                return 0;
        }


        if (npart <= 0) {
          PyErr_SetString( PyExc_ValueError, "Number of particles has to be larger than zero" );
          return 0;
        }


        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        dims[3] = nspecies+1;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 4, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny*nz;
        
        memset( grid, 0, cells*(nspecies+1)*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;
        cellsizez = bz / nz;
        
        mincellsize = min( cellsizex, min( cellsizey, cellsizez ) );

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        data_abund = (double*)PyArray_DATA(abund);
        
        xnuc = (double*)malloc( nspecies * sizeof(double) );
        lost = 0;
        mlost = 0;
        
        for (part=0; part<npart; part++) {
                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                if (h < 0.5*mincellsize) h = 0.5*mincellsize;
                h2 = h*h;

                m = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);
                
                for (i=0; i<nspecies; i++) {
                        xnuc[i] = *data_abund;
                        data_abund = (double*)((char*)data_abund + PyArray_STRIDES(abund)[1]);
                }
                data_abund = (double*)((char*)data_abund - nspecies*PyArray_STRIDES(abund)[1] + PyArray_STRIDES(abund)[0]);
                
                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                zmin = max( floor( (pz - h - cz + 0.5*bz) / cellsizez ), 0 );
                zmax = min( ceil( (pz + h - cz + 0.5*bz) / cellsizez ), nz-1 );

                zmid = floor( 0.5 * (zmin+zmax) + 0.5 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && zmin < nz && zmax >= 0) {
                        sum = 0.;
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                
                                        for (z0=zmid; z0>=zmin; z0--) {
                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                sum += 1.0;
                                        }

                                        for (z1=zmid+1; z1<=zmax; z1++) {
                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                sum += 1.0;
                                        }
                                }
                        }
                        
                        if (sum == 0.) {
                                lost++;
                                mlost += m;
                                continue;
                        }
                        
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                                                        
                                        for (z0=zmid; z0>=zmin; z0--) {
                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                grid[((x*ny + y)*nz + z0)*(nspecies+1)+nspecies] += m / sum;
                                                for(i=0; i<nspecies;i ++)
                                                        grid[((x*ny + y)*nz + z0)*(nspecies+1)+i] += xnuc[i] * m / sum;
                                        }

                                        for (z1=zmid+1; z1<=zmax; z1++) {
                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                grid[((x*ny + y)*nz + z1)*(nspecies+1)+nspecies] += m / sum;
                                                for(i=0; i<nspecies;i ++)
                                                        grid[((x*ny + y)*nz + z1)*(nspecies+1)+i] += xnuc[i] * m / sum;
                                        }
                                }
                        }       
                }
        }
        
        free( xnuc );
        
        mtot = 0;
        for(cell=0; cell<cells; cell++) {
                if(grid[cell*nspecies] > 0) {
                        m = grid[cell*(nspecies+1)];
                        for (i=0; i<nspecies; i++)
                                grid[cell*(nspecies+1)+i] /= m;
                        mtot += m;
                }
                
        }             

        printf( "Lost %d particles with a total mass of %g, mass on grid=%g.\n", lost, mlost, mtot );
        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _gatherAbundGrid(PyObject *self, PyObject *args, PyObject *kwargs) {
        PyArrayObject *pos, *mass, *rho, *abund, *pyGrid, *pyDensityField, *pyCoords;
        t_sph_tree tree;
        int npart, nx, ny, nz, nspecies, nneighbours, id;
        size_t cells;
        npy_intp dims[4];
        double bx, by, bz, cx, cy, cz;
        double *data_mass, *real_pos;
        double *grid;
        float *fgrid;
        double posz, poszmin, poszmax;
        long i, j;
        int gcheck, forceneighbourcount, is2d, converged, single_precision;
        time_t start;
        double *abundsum, masssum, *abundmin, *abundmax, *dens, maxratio, ratio, densitycut;
        int *specieslist, iratio, store, dummy, useCoords, nCoord, numthreads;
        FILE *fcheck;
        char *kwlist[] = {"pos", "mass", "rho", "abund", "nneighbours", "nx", "ny", "nz", "boxx", "boxy", "boxz", "centerx", "centery", "centerz", "forceneighbourcount", "densitycut", "densityfield", "gradientcheck", "single_precision", "coord", "numthreads", NULL};
        
        start = time( NULL );

        gcheck = 0;
        pyDensityField = NULL;
        densitycut = 0;
        forceneighbourcount = 0;
        single_precision = 0;
        pyCoords = NULL;
        numthreads = 1;
        if (!PyArg_ParseTupleAndKeywords( args, kwargs, "O!O!O!O!iiiidddddd|idO!iiO!i:gatherAbundGrid(pos, mass, rho, abund, nneighbours, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz, forceneighbourcount=0, densitycut=0, densityfield=None, gradientcheck=0, single_precision=False, coord=None, numthreads=1)", kwlist, &PyArray_Type, &pos, &PyArray_Type, &mass, &PyArray_Type, &rho, &PyArray_Type, &abund, &nneighbours, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz, &forceneighbourcount, &densitycut, &PyArray_Type, &pyDensityField, &gcheck, &single_precision, &PyArray_Type, &pyCoords, &numthreads )) {
                return 0;
        }

        printf( "This is gatherAbundGrid.\n" );

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE
            || PyArray_DESCR(pos)->byteorder != '=') {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3], type double and native byte order" );
                return 0;
        }
        npart = PyArray_DIMS(pos)[0];

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE || PyArray_DESCR(mass)->byteorder != '=') {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n], type double and native byte order" );
                return 0;
        }

        if (PyArray_NDIM(rho) != 1 || PyArray_TYPE(rho) != NPY_DOUBLE || PyArray_DESCR(rho)->byteorder != '=') {
                PyErr_SetString( PyExc_ValueError, "rho has to be of dimension [n], type double and native byte order" );
                return 0;
        }

        if (PyArray_NDIM(abund) != 2 || PyArray_TYPE(abund) != NPY_DOUBLE || PyArray_DESCR(abund)->byteorder != '=') {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n,nspecies], type double and native byte order" );
                return 0;
        }
        nspecies = PyArray_DIMS(abund)[1];
        
        if (npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(rho)[0] || npart != PyArray_DIMS(abund)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, rho and abund have to have the same size in the first dimension" );
                return 0;
        }

        printf( "%d particles with %d species found.\n", npart, nspecies );

        if (densitycut > 0.) {
          if (!pyDensityField ||
              PyArray_NDIM(pyDensityField) != 3 || 
              PyArray_DIMS(pyDensityField)[0] != nx || 
              PyArray_DIMS(pyDensityField)[1] != ny ||
              PyArray_DIMS(pyDensityField)[2] != nz ||
              PyArray_TYPE(pyDensityField) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "densityfield has to be of dimension [nx,ny,nz] and type double" );
                return 0;
          }

          printf( "Using densitycut: %g g/ccm\n", densitycut );
        }
        
        if(pyCoords) {
                if(PyArray_NDIM(pyCoords) != 2 || PyArray_DIMS(pyCoords)[1] != 3 || PyArray_TYPE(pyCoords) != NPY_DOUBLE) {
                        PyErr_SetString( PyExc_ValueError, "coord has to be of dimension [n,3] and type double" );
                        return 0;
                }
                nCoord = PyArray_DIMS(pyCoords)[0];
                useCoords = 1;
        }
        else {
                useCoords = 0;
        }
        
        poszmin =  1e200;
        poszmax = -1e200;
        for (id=0; id<npart; id++) {
                posz = *(double*)((char*)PyArray_DATA(pos) + id*PyArray_STRIDES(pos)[0] + 2*PyArray_STRIDES(pos)[1]);
                if (posz < poszmin) poszmin = posz;
                if (posz > poszmax) poszmax = posz;
        }
        
        if (poszmin == poszmax) {
                is2d = 1;
                printf( "Using 2d mode.\n" );
                if (!forceneighbourcount) {
                  printf( "Using exact number of neighbours!\n" );
                  forceneighbourcount = 1;
                }
        } else {
                is2d = 0;
        }
        
        
        specieslist = NULL;
        if (gcheck) {
                abundsum = (double*)malloc( nspecies * sizeof(double) );
                memset( abundsum, 0, nspecies * sizeof(double) );
                masssum = 0;

                for (id=0; id<npart; id++) {
                        double m = *(double*)((char*)PyArray_DATA(mass) + id*PyArray_STRIDES(mass)[0]);
                        masssum += m;
                        
                        int species;
                        for (species=0; species<nspecies; species++) {
                                abundsum[species] += *(double*)((char*)PyArray_DATA(abund) + id*PyArray_STRIDES(abund)[0] + species*PyArray_STRIDES(abund)[1]) * m;
                        }
                }

                int species;
                for (species=0; species<nspecies; species++) {
                        abundsum[species] /= masssum;
                }

                specieslist = (int*)malloc( gcheck * sizeof(int) );
                for (i=0; i<gcheck; i++) {
                        j = 0;
                        while ((j < i) && (abundsum[specieslist[j]] > abundsum[i])) j++;

                        store = i;
                        while (j <= i) {
                                dummy = specieslist[j];
                                specieslist[j] = store;
                                store = dummy;
                                j++;
                        }
                }

                for (i=gcheck; i<nspecies; i++) {
                        j = gcheck;

                        while ((j > 0) && (abundsum[specieslist[j-1]] < abundsum[i])) j--;
                        
                        if (j < gcheck) {
                                store = i;
                                while (j < gcheck) {
                                        dummy = specieslist[j];
                                        specieslist[j] = store;
                                        store = dummy;
                                        j++;
                                }
                        }
                }

                abundmin = (double*)malloc( gcheck * sizeof(double) );
                abundmax = (double*)malloc( gcheck * sizeof(double) );

                fcheck = fopen( "gradientcheck.dat", "w" );
                if(!fcheck) {
                        printf( "Failed to open file gradientcheck.dat.\n" );
                        return 0;
                }
                fwrite( &gcheck, sizeof(int), 1, fcheck );
                fwrite( specieslist, sizeof(int), gcheck, fcheck );

                free( abundsum );
        }

        if(useCoords) {
                dims[0] = nCoord;
                dims[1] = nspecies+1;
                printf( "Allocating %d grid for %d species, total size=%g GB.\n", (int)dims[0], (int)dims[1],
                        1e-9*dims[0]*dims[1]*8 );
                pyGrid = (PyArrayObject *)PyArray_SimpleNew( 2, dims, single_precision ? NPY_FLOAT : NPY_DOUBLE );
                if (!pyGrid) {
                        PyErr_SetString( PyExc_MemoryError, "Could not allocate memory for grid array." );
                        return 0;
                }
                
                nx = nCoord;
                ny = 1;
                nz = 1;
        } else {
                dims[0] = nx;
                dims[1] = ny;
                dims[2] = nz;
                dims[3] = nspecies+1;
                printf( "Allocating %dx%dx%d grid for %d species, total size=%g GB.\n", (int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3],
                        1e-9*dims[0]*dims[1]*dims[2]*dims[3]*8 );
                pyGrid = (PyArrayObject *)PyArray_SimpleNew( 4, dims, single_precision ? NPY_FLOAT : NPY_DOUBLE );
                if (!pyGrid) {
                        PyErr_SetString( PyExc_MemoryError, "Could not allocate memory for grid array." );
                        return 0;
                }
        }       

        cells = nx*ny*nz;
        if (single_precision) {
                fgrid = (float*)PyArray_DATA(pyGrid);
                memset( fgrid, 0, cells*(nspecies+1)*sizeof(*fgrid) );
        }
        else {
                grid = (double*)PyArray_DATA(pyGrid);
                memset( grid, 0, cells*(nspecies+1)*sizeof(*grid) );
        }

        
        if (forceneighbourcount) {
          dens = malloc( cells * sizeof(double) );
          memset( dens, 0, cells * sizeof(double) );
        }
        
        set_num_threads(numthreads);
        
        data_mass = (double*)PyArray_DATA(mass);
        
        real_pos = (double*)malloc( 3 * npart * sizeof( double ) );
        
        #pragma omp parallel for private(i, j)
        for (i=0; i<npart; i++) {
          for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pos) + i*PyArray_STRIDES(pos)[0] + j*PyArray_STRIDES(pos)[1]);
          }
        }

        createTree( &tree, npart, real_pos );

        size_t count_notconverged_all = 0;

        #pragma omp parallel
        {
        size_t x, y, z, cell;
        double cpx, cpy, cpz;
        int neighbourcount, part, nneighbours_real, species;
        double grid_pos[3];
        double h2, grid_hsml, grid_rho, r2, wk, r, a, m;
        
        int nextoutput = cells / 100;
        size_t count_notconverged = 0;
        int thread_id = get_thread_id();
        int slowly = 0;
        
        int *neighbours = malloc( tree.usednodes * sizeof(int) );
        
        r2 = 0;
        
        #pragma omp for schedule(dynamic) nowait reduction(+:count_notconverged_all)
        for (x=0; x<(size_t)nx; x++) {
                cpx = -0.5*bx + bx*(x+0.5)/nx + cx;
                for (y=0; y<(size_t)ny; y++) {
                        cpy = -0.5*by + by*(y+0.5)/ny + cy;
                        for (z=0; z<(size_t)nz; z++) {
                                cpz = -0.5*bz + bz*(z+0.5)/nz + cz;
                                
                                cell = (x*ny + y)*nz + z;

                                if (densitycut > 0 &&
                                    *(double*)((char*)PyArray_DATA(pyDensityField)  + x*PyArray_STRIDES(pyDensityField)[0] + 
                                               y*PyArray_STRIDES(pyDensityField)[1] + z*PyArray_STRIDES(pyDensityField)[2]) <= densitycut) {
                                      continue;
                                }
                                
                                if(useCoords) {
                                        grid_pos[0] = *(double*)((char*)PyArray_DATA(pyCoords) + x*PyArray_STRIDES(pyCoords)[0] + 0*PyArray_STRIDES(pyCoords)[1]);
                                        grid_pos[1] = *(double*)((char*)PyArray_DATA(pyCoords) + x*PyArray_STRIDES(pyCoords)[0] + 1*PyArray_STRIDES(pyCoords)[1]);
                                        grid_pos[2] = *(double*)((char*)PyArray_DATA(pyCoords) + x*PyArray_STRIDES(pyCoords)[0] + 2*PyArray_STRIDES(pyCoords)[1]);
                                } else {
                                        grid_pos[0] = cpx;
                                        grid_pos[1] = cpy;
                                        grid_pos[2] = cpz;
                                }
                                
                                if (forceneighbourcount) {
                                        grid_hsml = getNNeighbours( &tree, grid_pos, real_pos, nneighbours, &nneighbours_real, &neighbours, &converged );
                                        neighbourcount = nneighbours_real;

                                        if (!converged) count_notconverged++;
                                } else {
                                        grid_hsml = 0;
                                        calcHsml( &tree, grid_pos, real_pos, data_mass, nneighbours, &grid_hsml, &grid_rho );
                                        neighbourcount = getNeighbours( &tree, grid_pos, real_pos, grid_hsml, &neighbours );
                                }
                                h2 = grid_hsml * grid_hsml;

                                if (gcheck) {
                                        for (i=0; i<gcheck; i++) {
                                                abundmin[i] = 1.0;
                                                abundmax[i] = 0.0;
                                        }
                                }
                                
                                for (part=0; part<neighbourcount; part++) {
                                        int id = neighbours[part];
                                        
                                        r2 = (grid_pos[0]-real_pos[id*3+0])*(grid_pos[0]-real_pos[id*3+0])
                                           + (grid_pos[1]-real_pos[id*3+1])*(grid_pos[1]-real_pos[id*3+1])
                                           + (grid_pos[2]-real_pos[id*3+2])*(grid_pos[2]-real_pos[id*3+2]);
                                          
                                        if (r2 < h2) {
                                                wk = _getkernel( grid_hsml, r2 );
                                                
                                                /* if 2d => correct normalisation */
                                                if (is2d) wk *= grid_hsml;
                                                
                                                m = *(double*)((char*)PyArray_DATA(mass) + id*PyArray_STRIDES(mass)[0]);
                                                /*r = *(double*)((char*)PyArray_DATA(rho) + id*PyArray_STRIDES(rho)[0]);*/
                                                
#define tmp(g) {                                                        \
                g[((x*ny + y)*nz + z)*(nspecies+1)+nspecies] += wk * m; /* density */ \
                                                                        \
                for (species=0; species<nspecies; species++) {          \
                        a = *(double*)((char*)PyArray_DATA(abund) + id*PyArray_STRIDES(abund)[0] + species*PyArray_STRIDES(abund)[1]); \
                                                                        \
                        g[((x*ny + y)*nz + z)*(nspecies+1) + species] += wk * m * a; /* abundances */ \
                }                                                       \
                                                }
                                                if (single_precision) tmp(fgrid)
                                                else tmp(grid)
#undef tmp
                                                
                                                if (forceneighbourcount) {
                                                        r = *(double*)((char*)PyArray_DATA(rho) + id*PyArray_STRIDES(rho)[0]);
                                                        dens[ (x*ny + y)*nz + z ] += wk * m * r;
                                                }
                                                
                                                if (gcheck) {
                                                        for (i=0; i<gcheck; i++) {
                                                                species = specieslist[i];
                                                                a = *(double*)((char*)PyArray_DATA(abund) + id*PyArray_STRIDES(abund)[0] + species*PyArray_STRIDES(abund)[1]);
                                                                if (a < abundmin[i]) abundmin[i] = a;
                                                                if (a > abundmax[i]) abundmax[i] = a;
                                                        }
                                                }
                                        }                                       
                                }
                                
                                if (gcheck && neighbourcount) {
                                        maxratio = 0;
                                        iratio = -1;
                                        for (i=0; i<gcheck; i++) {
                                                if (abundmax[i] > 1e-3) {
                                                        ratio = (abundmax[i]-abundmin[i])/(abundmax[i]+abundmin[i]);
                                                        if ( ratio > maxratio ) {
                                                                maxratio = ratio;
                                                                iratio = i;
                                                        }
                                                }
                                        }

                                        if (maxratio > 0.1) {
                                                printf( "r: %g km, massratio: %g, min: %g, max: %g\n", sqrt(r2) / 1e5, maxratio, abundmin[iratio], abundmax[iratio] );
                                                fwrite( &maxratio, sizeof(double), 1, fcheck );
                                                fwrite( &grid_pos, sizeof(double), 3, fcheck );
                                                fwrite( abundmin, sizeof(double), gcheck, fcheck );
                                                fwrite( abundmax, sizeof(double), gcheck, fcheck );

                                                fwrite( &neighbourcount, sizeof(int), 1, fcheck );
                                                fwrite( neighbours, sizeof(int), neighbourcount, fcheck );
                                        }
                                }         
                                
                                if( thread_id == numthreads - 1 ) {
                                        if (cell >= (size_t)nextoutput) {
                                          time_t now = time( NULL );
                                          double runtime = difftime( now, start );
                                  
                                          if ((size_t)nextoutput == cells / 100) {
                                            if ( runtime > 60. ) {
                                              slowly = 1;
                                            } else {
                                              nextoutput = cells / 10;
                                            }
                                          }

                                          printf( "%zd / %zd cells done (%d%%): %ds elapsed, ~%ds remaining\n", cell, cells, (int)floor(100.0*(double)cell/(double)cells), (int)(runtime), (int)(runtime/cell*(cells-cell)) );

                                          if (slowly)
                                            nextoutput += cells / 100;
                                          else
                                            nextoutput += cells /  10;
                                        }
                                }
                        }
                }
        }

        count_notconverged_all += count_notconverged;

        free(neighbours);
        } // close parallel section
        
        if (count_notconverged_all) {
                printf( "Neighbour search not converged for %zd cells => be careful\n", count_notconverged_all );
        }

        if (gcheck) {
                fclose( fcheck );
                free( specieslist );
                free( abundmin );
                free( abundmax );
        }
        
        free( real_pos );
        
        {
        size_t cell;
        int species;
        double r;
                
#define tmp(g)                                                          \
        for (cell=0; cell<cells; cell++) {                              \
                r = g[cell*(nspecies+1)+nspecies];                      \
                                                                        \
                if (r > 0) {                                            \
                        for (species=0; species<nspecies; species++) {  \
                                g[cell*(nspecies+1)+species] /= r;      \
                        }                                               \
                                                                        \
                        if (forceneighbourcount) g[cell*(nspecies+1)+nspecies] = dens[cell] / r; \
                }                                                       \
        }
        if (single_precision) tmp(fgrid)
        else tmp(grid)
#undef tmp
        }

        if (forceneighbourcount) free( dens );

        freeTree( &tree );

        time_t now = time( NULL );
        double runtime = difftime( now, start );
        printf( "Calculation took %ds\n", (int)runtime );
        return PyArray_Return( pyGrid );
}

PyObject* _calcSlice(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *rho, *value, *pyGrid;
        int npart, nx, ny, cells;
        int dims[2];
        double bx, by, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass, *data_rho, *data_value;
        double *grid;
        int part;
        double px, py, pz, h, m, r, v, cpx, cpy, r2, h2;
        double p[3];
        int x, y;
        int xmin, xmax, ymin, ymax, axis0, axis1;
        double cellsizex, cellsizey;
        clock_t start;
        
        start = clock();

        axis0 = 0;
        axis1 = 1;
        if (!PyArg_ParseTuple( args, "O!O!O!O!O!iiddddd|ii:calcSlice( pos, hsml, mass, rho, value, nx, ny, boxx, boxy, centerx, centery, centerz, [axis0, axis1] )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &rho, &PyArray_Type, &value, &nx, &ny, &bx, &by, &cx, &cy, &cz, &axis0, &axis1 )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(rho) != 1 || PyArray_TYPE(rho) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "rho has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(rho)[0] || npart != PyArray_DIMS(value)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml, mass, rho and value have to have the same size in the first dimension" );
                return 0;
        }
        dims[0] = nx;
        dims[1] = ny;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny;
        memset( grid, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        data_rho = (double*)PyArray_DATA(rho);
        data_value = (double*)PyArray_DATA(value);

        for (part=0; part<npart; part++) {
                p[0] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[1] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[2] = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                px = p[ axis0 ];
                py = p[ axis1 ];
                pz = p[ 3 - axis0 - axis1 ];
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                m = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                r = *data_rho;
                data_rho = (double*)((char*)data_rho + PyArray_STRIDES(rho)[0]);

                v = *data_value;
                data_value = (double*)((char*)data_value + PyArray_STRIDES(value)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && abs(pz-cz) < h) {
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        r2 = sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cz);
                                        if (r2 > h2) continue;
                                        grid[x*ny + y] += _getkernel( h, r2 ) * m * v / r;
                                }
                        }       
                }
        }

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _calcGridMassWeight(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *value, *pyGridMass, *pyGridValue;
        int npart, nx, ny, nz, cells;
        int dims[3];
        double bx, by, bz, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass, *data_value;
        double *gridmass, *gridvalue, *massend, *massiter, *valueiter;
        int part;
        double px, py, pz, h, h2, m, v, cpx, cpy, cpz, r2, dmass;
        int x, y, z0, z1;
        int xmin, xmax, ymin, ymax, zmin, zmax, zmid;
        double cellsizex, cellsizey, cellsizez;

        if (!PyArg_ParseTuple( args, "O!O!O!O!O!iiidddddd:calcGridMassWeight( pos, hsml, mass, value, massgrid, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &value, &PyArray_Type, &pyGridMass, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0] || npart != PyArray_DIMS(value)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and value have to have the same size in the first dimension" );
                return 0;
        }

        if (PyArray_NDIM(pyGridMass) != 3 || PyArray_DIMS(pyGridMass)[0] != nx || PyArray_DIMS(pyGridMass)[1] != ny || PyArray_DIMS(pyGridMass)[2] != nz) {
                PyErr_SetString( PyExc_ValueError, "massgrid has to have 3 dimensions: [nx,ny,nz]" );
                return 0;
        }

        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        cells = nx*ny*nz;
        
        gridmass = (double*)PyArray_DATA(pyGridMass);
        
        pyGridValue = (PyArrayObject *)PyArray_FromDims( 3, dims, NPY_DOUBLE );
        gridvalue = (double*)PyArray_DATA(pyGridValue);
        memset( gridvalue, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;
        cellsizez = bz / nz;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        data_value = (double*)PyArray_DATA(value);

        for (part=0; part<npart; part++) {
                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                m = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                v = *data_value;
                data_value = (double*)((char*)data_value + PyArray_STRIDES(value)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                zmin = max( floor( (pz - h - cz + 0.5*bz) / cellsizez ), 0 );
                zmax = min( ceil( (pz + h - cz + 0.5*bz) / cellsizez ), nz-1 );

                zmid = floor( 0.5 * (zmin+zmax) + 0.5 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && zmin < nz && zmax >= 0) {
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        for (z0=zmid; z0>=zmin; z0--) {
                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                dmass = _getkernel( h, r2 ) * m;
                                                gridvalue[(x*ny + y)*nz + z0] += dmass * v;
                                        }

                                        for (z1=zmid+1; z1<=zmax; z1++) {
                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                dmass = _getkernel( h, r2 ) * m;
                                                gridvalue[(x*ny + y)*nz + z1] += dmass * v;
                                        }
                                }
                        }       
                }
        }
        
        massend = &gridmass[ cells ];
        for (massiter = gridmass, valueiter = gridvalue; massiter != massend; massiter++, valueiter++) {
                if (*massiter > 0)
                        *valueiter /= *massiter;
        }
        
        return PyArray_Return( pyGridValue );
}

PyObject* _calcDensGrid(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *pyGrid;
        int npart, nx, ny, nz, cells;
        int dims[3];
        double bx, by, bz, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass;
        double *grid;
        int part;
        double px, py, pz, h, h2, v, cpx, cpy, cpz, r2;
        int x, y, z0, z1;
        int xmin, xmax, ymin, ymax, zmin, zmax, zmid;
        double cellsizex, cellsizey, cellsizez;
        clock_t start;
        
        start = clock();
        
        if (!PyArg_ParseTuple( args, "O!O!O!iiidddddd:calcDensGrid( pos, hsml, mass, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and mass have to have the same size in the first dimension" );
                return 0;
        }

        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 3, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny*nz;
        memset( grid, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;
        cellsizez = bz / nz;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        
        for (part=0; part<npart; part++) {
                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                v = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                zmin = max( floor( (pz - h - cz + 0.5*bz) / cellsizez ), 0 );
                zmax = min( ceil( (pz + h - cz + 0.5*bz) / cellsizez ), nz-1 );

                zmid = floor( 0.5 * (zmin+zmax) + 0.5 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && zmin < nz && zmax >= 0) {
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        for (z0=zmid; z0>=zmin; z0--) {
                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                grid[(x*ny + y)*nz + z0] += _getkernel( h, r2 ) * v;
                                        }

                                        for (z1=zmid+1; z1<=zmax; z1++) {
                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                grid[(x*ny + y)*nz + z1] += _getkernel( h, r2 ) * v;
                                        }
                                }
                        }       
                }
        }

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _calcDensSlice(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *pyGrid;
        int npart, nx, ny, cells;
        int dims[2];
        double bx, by, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass;
        double *grid;
        int part;
        double px, py, pz, h, v, cpx, cpy, r2, h2;
        double p[3];
        int x, y;
        int xmin, xmax, ymin, ymax, axis0, axis1;
        double cellsizex, cellsizey;
        clock_t start;
        
        start = clock();

        axis0 = 0;
        axis1 = 1;
        if (!PyArg_ParseTuple( args, "O!O!O!iiddddd|ii:calcDensSlice( pos, hsml, mass, nx, ny, boxx, boxy, centerx, centery, centerz, [axis0, axis1] )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &nx, &ny, &bx, &by, &cx, &cy, &cz, &axis0, &axis1 )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and mass have to have the same size in the first dimension" );
                return 0;
        }

        dims[0] = nx;
        dims[1] = ny;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny;
        memset( grid, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);

        for (part=0; part<npart; part++) {
                p[0] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[1] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[2] = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                px = p[ axis0 ];
                py = p[ axis1 ];
                pz = p[ 3 - axis0 - axis1 ];
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                v = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && abs(pz-cz) < h) {
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        r2 = sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cz);
                                        if (r2 > h2) continue;
                                        grid[x*ny + y] += _getkernel( h, r2 ) * v;
                                }
                        }       
                }
        }

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _calcCylinderAverage(PyObject *self, PyObject *args) {
        PyArrayObject *pyGrid, *pyNewgrid;
        int dims[2], cells;
        double *newgrid, *count;
        int x, y, z, nx, ny, nz, nr, r;

        if (!PyArg_ParseTuple( args, "O!:calcCylinderAverage( grid )", &PyArray_Type, &pyGrid )) {
                return 0;
        }

        if (PyArray_NDIM(pyGrid) != 3 || PyArray_TYPE(pyGrid) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "grid has to be of dimensions [nx,ny,nz] and type double" );
                return 0;
        }

        nx = PyArray_DIMS(pyGrid)[0];
        ny = PyArray_DIMS(pyGrid)[1];
        nz = PyArray_DIMS(pyGrid)[2];
        nr = min( ny, nz );
        
        dims[0] = nx;
        dims[1] = nr;
        cells = nx*nr;
        pyNewgrid = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
        newgrid = (double*)PyArray_DATA(pyNewgrid);
        memset( newgrid, 0, cells*sizeof(double) );

        count = (double*)malloc( cells*sizeof(double) );
        memset( count, 0, cells*sizeof(double) );

        for (x=0; x<nx; x++) for (y=0; y<ny; y++) for (z=0; z<nz; z++) {
                r = floor( sqrt( sqr(y-ny/2.0+0.5) + sqr(z-nz/2.0+0.5) ) );
                if (r >= nr/2) continue;

                newgrid[x*nr+r+nr/2] += *(double*)( PyArray_DATA(pyGrid) + PyArray_STRIDES(pyGrid)[0]*x + PyArray_STRIDES(pyGrid)[1]*y + PyArray_STRIDES(pyGrid)[2]*z );
                count[x*nr+r+nr/2] += 1;
        }
        
        for (x=0; x<nx; x++) for (r=0; r<nr/2; r++) {
                if (count[x*nr+r+nr/2] > 0) {
                        newgrid[x*nr+r+nr/2] /= count[x*nr+r+nr/2];
                        newgrid[x*nr-r+nr/2-1] = newgrid[x*nr+r+nr/2];
                }
        }

        free( count );
        return PyArray_Return( pyNewgrid );
}

PyObject* _calcRadialProfile(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *data, *pyProfile;
        int npart, nshells, mode;
        int dims[2];
        int *count;
        double cx, cy, cz, dr;
        double *data_pos, *data_data;
        double *profile;
        int part, shell;
        double px, py, pz, d, rr, v;
        clock_t start;
        
        start = clock();

        mode = 1;
        nshells = 200;
        dr = 0;
        cx = cy = cz = 0;
        if (!PyArg_ParseTuple( args, "O!O!|iidddd:calcRadialProfile( pos, data, mode, nshells, dr, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &data, &mode, &nshells, &dr, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(data) != 1 || PyArray_TYPE(data) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "data has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(data)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos and data have to have the same size in the first dimension" );
                return 0;
        }
        dims[0] = 2;
        dims[1] = nshells;
        pyProfile = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
        profile = (double*)PyArray_DATA(pyProfile);
        memset( profile, 0, 2*nshells*sizeof(double) );

        count = (int*)malloc( nshells*sizeof(int) );
        memset( count, 0, nshells*sizeof(int) );

        if (!dr) {
                data_pos = (double*)PyArray_DATA(pos);
                for (part=0; part<npart; part++) {
                        px = *data_pos;
                        data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                        py = *data_pos;
                        data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                        pz = *data_pos;
                        data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);

                        rr = sqrt( sqr(px-cx) + sqr(py-cy) + sqr(pz-cz) );
                        if (rr > dr)
                                dr = rr;
                }
                dr /= nshells;
        }

        data_pos = (double*)PyArray_DATA(pos);
        data_data = (double*)PyArray_DATA(data);

        for (part=0; part<npart; part++) {

                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);

                d = *data_data;
                data_data = (double*)((char*)data_data + PyArray_STRIDES(data)[0]);

                rr = sqrt( sqr(px-cx) + sqr(py-cy) + sqr(pz-cz) );
                shell = floor( rr / dr );

                if (shell < nshells) {
                        profile[ shell ] += d;
                        count[ shell ] += 1;
                }
        }

        for (shell=0; shell<nshells; shell++) {
                profile[ nshells + shell ] = dr * (shell + 0.5);
        }

        switch (mode) {
                // sum
                case 0:
                        break;
                // density
                case 1:
                        for (shell=0; shell<nshells; shell++) {
                                v = 4.0 / 3.0 * M_PI * dr*dr*dr * ( ((double)shell+1.)*((double)shell+1.)*((double)shell+1.) - (double)shell*(double)shell*(double)shell );
                                profile[shell] /= v;
                        }
                        break;
                // average
                case 2:
                        for (shell=0; shell<nshells; shell++) if (count[shell] > 0) profile[shell] /= count[shell];
                        break;
        }

        free( count );

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyProfile );
}

PyObject* _calcAbundGridSph(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *abund, *pyGrid;
        int npart, nx, ny, nz, nspecies;
        size_t cell, cells;
        int dims[4];
        double bx, by, bz, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass, *data_abund;
        double *grid;
        int part, species;
        double *xnuc;
        double px, py, pz, h, h2, m, cpx, cpy, cpz, r2, kk, rho;
        size_t x, y, z0, z1;
        int xmin, xmax, ymin, ymax, zmin, zmax, zmid;
        int slowly, nextoutput;
        double cellsizex, cellsizey, cellsizez, runtime;
        clock_t start;
        
        start = clock();

        if (!PyArg_ParseTuple( args, "O!O!O!O!iiidddddd:calcAbundGrid( pos, hsml, mass, abund, nx, ny, nz, boxx, boxy, boxz, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &abund, &nx, &ny, &nz, &bx, &by, &bz, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(abund) != 2 || PyArray_TYPE(abund) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "abund has to be of dimension [n,nspecies] and type double" );
                return 0;
        }
        
        nspecies = PyArray_DIMS(abund)[1];

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(abund)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and abund have to have the same size in the first dimension" );
                return 0;
        }
        
        xnuc = (double*)malloc( nspecies * sizeof( double ) );

        dims[0] = nx;
        dims[1] = ny;
        dims[2] = nz;
        dims[3] = nspecies+1;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 4, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny*nz;
        memset( grid, 0, cells*(nspecies+1)*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;
        cellsizez = bz / nz;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        data_abund = (double*)PyArray_DATA(abund);

        slowly = 0;
        nextoutput = npart / 100;

        for (part=0; part<npart; part++) {
                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                m = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                for (species = 0; species < nspecies; species++ ) {
                        xnuc[species] = *data_abund;
                        data_abund = (double*)((char*)data_abund + PyArray_STRIDES(abund)[1]);
                }
                data_abund = (double*)((char*)data_abund - nspecies*PyArray_STRIDES(abund)[1] + PyArray_STRIDES(abund)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                zmin = max( floor( (pz - h - cz + 0.5*bz) / cellsizez ), 0 );
                zmax = min( ceil( (pz + h - cz + 0.5*bz) / cellsizez ), nz-1 );

                zmid = floor( 0.5 * (zmin+zmax) + 0.5 );

                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && zmin < nz && zmax >= 0) {
                        for (x=xmin; x<=(size_t)xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=(size_t)ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        for (z0=zmid; z0>=(size_t)zmin; z0--) {
                                                cpz = -0.5*bz + bz*(z0+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                kk = _getkernel( h, r2 ) * m;
                                                for (species = 0; species < nspecies; species++ )
                                                        grid[((x*ny + y)*nz + z0)*(nspecies+1) + species] += kk * xnuc[species];
                                                grid[((x*ny + y)*nz + z0)*(nspecies+1) + nspecies] += kk;
                                        }

                                        for (z1=zmid+1; z1<=(size_t)zmax; z1++) {
                                                cpz = -0.5*bz + bz*(z1+0.5)/nz;
                                                r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                                                if (r2 > h2) break;
                                                
                                                kk = _getkernel( h, r2 ) * m;
                                                for (species = 0; species < nspecies; species++ )
                                                        grid[((x*ny + y)*nz + z1)*(nspecies+1) + species] += kk * xnuc[species];
                                                grid[((x*ny + y)*nz + z1)*(nspecies+1) + nspecies] += kk;
                                        }
                                }
                        }       
                }

                if (part >= nextoutput) {
                        runtime = ((double)clock()-(double)start)/CLOCKS_PER_SEC;
                                  
                        if (nextoutput == npart / 100) {
                          if ( runtime > 60. ) {
                            slowly = 1;
                          } else {
                            nextoutput = 0;
                          }
                        }

                        printf( "%d / %d particles done (%d%%): %ds elapsed, ~%ds remaining\n", part, npart, (int)floor(100.0*(double)part/(double)npart), (int)(runtime), (int)(runtime/part*(npart-part)) );

                        if (slowly)
                          nextoutput += npart / 100;
                        else
                          nextoutput += npart /  10;
                }
        }
        
        free( xnuc );

        for (cell=0; cell<cells; cell++) {              
                rho = grid[cell*(nspecies+1)+nspecies];
                if (rho > 0) {
                        for (species=0; species<nspecies; species++) {
                          grid[cell*(nspecies+1)+species] /= rho;
                        }                                               
                }                                                       
        }

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _calcAbundSphere(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *abund, *pyGrid;
        int npart, nradius, ntheta, nphi, cells, nspecies;
        int dims[4];
        double radius, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass, *data_abund;
        double *grid;
        int part, species;
        double *xnuc;
        double px, py, pz, h, h2, m, r, dr, cpx, cpy, cpz, r2, kk;
        double vr, vtheta, vphi;
        int ir, itheta, iphi;
        int minradius, maxradius;
        clock_t start;
        
        start = clock();

        if (!PyArg_ParseTuple( args, "O!O!O!O!iiidddd:calcAbundSphere( pos, hsml, mass, abund, nradius, ntheta, nphi, radius, centerx, centery, centerz )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &PyArray_Type, &abund, &nradius, &ntheta, &nphi, &radius, &cx, &cy, &cz )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(abund) != 2 || PyArray_TYPE(abund) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "abund has to be of dimension [n,nspecies] and type double" );
                return 0;
        }
        
        nspecies = PyArray_DIMS(abund)[1];

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]  || npart != PyArray_DIMS(abund)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and abund have to have the same size in the first dimension" );
                return 0;
        }
        
        xnuc = (double*)malloc( nspecies * sizeof( double ) );

        dims[0] = nradius;
        dims[1] = ntheta;
        dims[2] = nphi;
        dims[3] = nspecies+1;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 4, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nradius*ntheta*nphi*(nspecies+1);
        memset( grid, 0, cells*sizeof(double) );

        dr = radius / nradius;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);
        data_abund = (double*)PyArray_DATA(abund);
                
        for (part=0; part<npart; part++) {
                px = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                py = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                pz = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                m = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                for (species = 0; species < nspecies; species++ ) {
                        xnuc[species] = *data_abund;
                        data_abund = (double*)((char*)data_abund + PyArray_STRIDES(abund)[1]);
                }
                data_abund = (double*)((char*)data_abund - nspecies*PyArray_STRIDES(abund)[1] + PyArray_STRIDES(abund)[0]);

                r = sqrt( px*px + py*py + pz*pz );
                minradius = max( 0, floor( (r-h-0.5) / dr ) );
                maxradius = min( nradius-1, floor( (r+h+0.5) / dr ) );

                for (ir=minradius; ir<=maxradius; ir++)
                for (itheta=0; itheta<ntheta; itheta++)
                for (iphi=0; iphi<nphi; iphi++) {
                        vr = radius * (ir+0.5) / nradius;
                        vtheta = M_PI * (itheta+0.5) / ntheta;
                        vphi = 2. * M_PI * (iphi+0.5) / nphi;
                        cpx = vr * sin( vtheta ) * cos( vphi );
                        cpy = vr * sin( vtheta ) * sin( vphi );
                        cpz = vr * cos( vtheta );

                        r2 = ( sqr(px-cpx-cx) + sqr(py-cpy-cy) + sqr(pz-cpz-cz) );
                        if (r2 > h2) continue;
                                                
                        kk = _getkernel( h, r2 ) * m;
                        for (species = 0; species < nspecies; species++ ) {
                                grid[((ir*ntheta + itheta)*nphi + iphi)*(nspecies+1) + species] += kk * xnuc[species];
                        }
                        grid[((ir*ntheta + itheta)*nphi + iphi)*(nspecies+1) + nspecies] += kk;
                }
        }
        
        free( xnuc );

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

PyObject* _calcASlice(PyObject *self, PyObject *args, PyObject *kwargs) {
        PyArrayObject *pos, *value, *grad, *pyGrid, *pyNeighbours, *pyContours, *pcenter;
        int npart, nx, ny, nz, axis0, axis1, proj, grid3D, ngbs, numthreads;
        int x, y, z, i, j, cell;                                                                               
        double bx, by, bz, cx, cy, cz;
        double *real_pos, *real_value, *real_grad, *real_pcenter;
        double *grid;
        double cellsizex, cellsizey, cellsizez;
        int neighbour, *neighbours, *contours;
        t_sph_tree tree;
        PyObject *dict;
        double start, starttree, end;
        char *kwlist[] = {"pos", "value", "nx", "ny", "boxx", "boxy", "centerx", "centery", "centerz", "axis0", "axis1", "proj", "grad", "pcenter", "nz", "boxz", "grid3D", "ngbs", "numthreads", NULL} ;
        
        start = get_time();

        axis0 = 0;
        axis1 = 1;
        proj = 0;
        grad = NULL;
        pcenter = NULL;
        pyNeighbours = NULL;
        pyContours = NULL;
        nz = 0;
        bz = 0;
        grid3D = 0;
        ngbs = 1;
        numthreads = 1;
        if (!PyArg_ParseTupleAndKeywords( args, kwargs, "O!O!iiddddd|iiiO!O!idiii:calcASlice( pos, value, nx, ny, boxx, boxy, centerx, centery, centerz, [axis0, axis1, proj, grad, pcenter, nz, boxz, grid3D, ngbs, numthreads] )", kwlist, &PyArray_Type, &pos, &PyArray_Type, &value, &nx, &ny, &bx, &by, &cx, &cy, &cz, &axis0, &axis1, &proj, &PyArray_Type, &grad, &PyArray_Type, &pcenter, &nz, &bz, &grid3D, &ngbs, &numthreads )) {
                return 0;
        }
        
        if (proj || grid3D)
        {
                if (nz == 0) {
                        nz = max( nx, ny );
                        printf( "nz=0, setting it to nz=%d.\n", nz );
                }
                if (bz == 0) {
                        bz = max( bx, by );
                        printf( "boxz=0, setting it to boxz=%g.\n", bz );
                }
        }
        else
        {
                nz = 1;
                bz = 0;
        }
        
        if (proj || grid3D)
                ngbs = 0;

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n] and type double" );
                return 0;
        }

        
        if (grad && (PyArray_NDIM(grad) != 2 || PyArray_DIMS(grad)[1] != 3 || PyArray_TYPE(grad) != NPY_DOUBLE)) {
                PyErr_SetString( PyExc_ValueError, "grad has to be of dimension [n,3] and type double" );
                return 0;
        }

        if (pcenter && (PyArray_NDIM(pcenter) != 2 || PyArray_DIMS(pcenter)[1] != 3 || PyArray_TYPE(pcenter) != NPY_DOUBLE)) {
          PyErr_SetString( PyExc_ValueError, "pcenter has to be of dimension [n,3] and type double" );
          return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(value)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos and value have to have the same size in the first dimension" );
                return 0;
        }
        
        if (grid3D) {
                npy_intp dims[3];
                dims[0] = nx;
                dims[1] = ny;
                dims[2] = nz;
                pyGrid = (PyArrayObject *)PyArray_SimpleNew( 3, (npy_intp*)dims, NPY_DOUBLE );
                grid = (double*)PyArray_DATA(pyGrid);
                memset( grid, 0, nx * ny * nz * sizeof(double) );
                printf( "Doing 3D Grid of size %d x %d x %d\n", nx, ny, nz );
        } else {
                npy_intp dims[2];
                dims[0] = nx;
                dims[1] = ny;
                pyGrid = (PyArrayObject *)PyArray_SimpleNew( 2, (npy_intp*)dims, NPY_DOUBLE );
                grid = (double*)PyArray_DATA(pyGrid);
                memset( grid, 0, nx * ny * sizeof(double) );
                
                if (ngbs) {
                        pyNeighbours = (PyArrayObject *)PyArray_SimpleNew( 2, (npy_intp*)dims, NPY_INT );
                        neighbours = (int*)PyArray_DATA(pyNeighbours);
                        pyContours = (PyArrayObject *)PyArray_SimpleNew( 2, (npy_intp*)dims, NPY_INT );
                        contours = (int*)PyArray_DATA(pyContours);
                        memset( contours, 0, nx * ny * sizeof(int) );
                }
        }

        set_num_threads(numthreads);

        real_pos = (double*)malloc( 3 * npart * sizeof( double ) );
        real_value = (double*)malloc( npart * sizeof( double ) );
        if (grad) real_grad = (double*)malloc( 3 * npart * sizeof( double ) );
        if (pcenter) real_pcenter = (double*)malloc( 3 * npart * sizeof( double ) );

        #pragma omp parallel for private(i, j)
        for (i=0; i<npart; i++) {
          for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pos) + i*PyArray_STRIDES(pos)[0] + j*PyArray_STRIDES(pos)[1]);
            if (grad) real_grad[i*3+j] = *(double*)((char*)PyArray_DATA(grad) + i*PyArray_STRIDES(grad)[0] + j*PyArray_STRIDES(grad)[1]);
            if (pcenter) real_pcenter[i*3+j] = *(double*)((char*)PyArray_DATA(pcenter) + i*PyArray_STRIDES(pcenter)[0] + j*PyArray_STRIDES(pcenter)[1]);
          }
          real_value[i] = *(double*)((char*)PyArray_DATA(value) + i*PyArray_STRIDES(value)[0]);
        }

        createTree( &tree, npart, real_pos );

        cellsizex = bx / nx;
        cellsizey = by / ny;

        if (proj || grid3D)
                cellsizez = bz / nz;
        else
                cellsizez = 0;

        neighbour = -1;
        cell = 0;

        printf("Starting tree walk with %d thread(s)\n", numthreads);

        starttree = get_time();

        #pragma omp parallel private(x, y, z) firstprivate(cell, neighbour)
        {
          int *worklist = malloc( tree.usednodes * sizeof(int) );

          double coord[3];
          int cnt = 0;
          int thread_id = get_thread_id();

          coord[ 3 - axis0 - axis1 ] = cz;

          #pragma omp for schedule(dynamic, 1) nowait
          for (x=0; x<nx; x++) {
            coord[ axis0 ] = cx - 0.5 * bx + cellsizex * (0.5 + x);
            for (y=0; y<ny; y++) {
              if(!grid3D) cell = y + ny * x;
              coord[ axis1 ] = cy - 0.5 * by + cellsizey * (0.5 + y);
              for (z=0; z<nz; z++) {
                if(grid3D) cell = z + nz * (y + ny * x);
                coord[ 3-axis0-axis1 ] = cz - 0.5 * bz + cellsizez * (0.5 + z);

                getNearestNeighbour( &tree, real_pos, coord, &neighbour, worklist );

                if (!grad)
                  {
                    #pragma omp atomic
                    grid[ cell ] += real_value[ neighbour ];
                  }
                else {
                  if(pcenter) {
                    #pragma omp atomic
                    grid[ cell ] += real_value[ neighbour ]
                      + (coord[0] - real_pcenter[ neighbour*3 + 0 ]) * real_grad[ neighbour*3 + 0 ]
                      + (coord[1] - real_pcenter[ neighbour*3 + 1 ]) * real_grad[ neighbour*3 + 1 ]
                      + (coord[2] - real_pcenter[ neighbour*3 + 2 ]) * real_grad[ neighbour*3 + 2 ];
                  } else {
                    #pragma omp atomic
                    grid[ cell ] += real_value[ neighbour ]
                      + (coord[0] - real_pos[ neighbour*3 + 0 ]) * real_grad[ neighbour*3 + 0 ]
                      + (coord[1] - real_pos[ neighbour*3 + 1 ]) * real_grad[ neighbour*3 + 1 ]
                      + (coord[2] - real_pos[ neighbour*3 + 2 ]) * real_grad[ neighbour*3 + 2 ];
                  }
                }
                if (ngbs) neighbours[ cell ] = neighbour;
              }
            }
            if (ngbs) neighbour = neighbours[ cell - ny ];

            if(thread_id == numthreads - 1)
              if(x / (nx/10) > cnt)
                {
                  cnt++;
                  printf("Done iter %4d of %4d\n", x+1, nx);
                }
          }
          
          free(worklist);
        }

        end = get_time();

        printf( "Tree walk took %gs\n", end-starttree );

        freeTree( &tree );

        free( real_pos );
        free( real_value );
        if (grad) free( real_grad );
        if (pcenter) free( real_pcenter );

        if (ngbs){
          #pragma omp parallel for private(x, y, neighbour)
          for (x=1; x<nx-1; x++) {
            for (y=1; y<ny-1; y++) {
              neighbour = neighbours[ x*ny + y ];
              if (neighbours[ (x-1)*ny + y-1 ] != neighbour ||
                  neighbours[  x   *ny + y-1 ] != neighbour ||
                  neighbours[ (x+1)*ny + y-1 ] != neighbour ||
                  neighbours[ (x-1)*ny + y   ] != neighbour ||
                  neighbours[ (x+1)*ny + y   ] != neighbour ||
                  neighbours[ (x-1)*ny + y+1 ] != neighbour ||
                  neighbours[  x   *ny + y+1 ] != neighbour ||
                  neighbours[ (x+1)*ny + y+1 ] != neighbour) {
                contours[ x*ny + y ] = 1;
              } else {
                contours[ x*ny + y ] = 0;
              }
            }
          }
        }

        dict = PyDict_New();
        PyDict_SetStolenItem( dict, "grid", (PyObject*)pyGrid );
        
        if (ngbs)
        {
                PyDict_SetStolenItem( dict, "neighbours", (PyObject*)pyNeighbours );
                PyDict_SetStolenItem( dict, "contours", (PyObject*)pyContours );
        }

        end = get_time();

        printf( "Calculation took %gs\n", end-start );

        return dict;
}

PyObject* _calcASlice2(PyObject *self, PyObject *args, PyObject *kwargs) {
        PyArrayObject *pos, *value, *pyGrid;
        int npart, nx, ny, axis0, axis1, proj, cic, tsc;
        int i, j, n;                                                                               
        double bx, by, cx, cy;
        double *real_pos, *real_value;
        double *grid;
        double cellsizex, cellsizey;
        PyObject *dict;
        clock_t start;
        char *kwlist[] = {"pos", "value", "nx", "ny", "boxx", "boxy", "centerx", "centery", "axis0", "axis1", "proj", "cic", "tsc", NULL} ;
        
        start = clock();

        axis0 = 0;
        axis1 = 1;
        proj = 0;
        if (!PyArg_ParseTupleAndKeywords( args, kwargs, "O!O!iidddd|iiiii:calcASlice2( pos, value, nx, ny, boxx, boxy, centerx, centery, [axis0, axis1, proj, cic, tsc] )", kwlist, &PyArray_Type, &pos, &PyArray_Type, &value, &nx, &ny, &bx, &by, &cx, &cy, &axis0, &axis1, &proj, &cic, &tsc )) {
                return 0;
        }
        
        if(!proj)
        {
                PyErr_SetString( PyExc_ValueError, "Only projection mode is allowed!" );
                return 0;
        }

        if(cic && tsc)
        {
                PyErr_SetString( PyExc_ValueError, "Please enable only one scheme between cic and tsc. Not enabling both is a valid option but results in nearest point assignment scheme" );
                return 0;
        }
        
        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(value) != 1 || PyArray_TYPE(value) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "value has to be of dimension [n] and type double" );
                return 0;
        }

        
        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(value)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos and value have to have the same size in the first dimension" );
                return 0;
        }
        
        if(proj) {
                npy_intp dims[2];
                dims[0] = nx;
                dims[1] = ny;
                pyGrid = (PyArrayObject *)PyArray_SimpleNew( 2, (npy_intp*)dims, NPY_DOUBLE );
                grid = (double*)PyArray_DATA(pyGrid);
                memset( grid, 0, nx * ny * sizeof(double) );
        }

        real_pos = (double*)malloc( 3 * npart * sizeof( double ) );
        real_value = (double*)malloc( npart * sizeof( double ) );
        for (i=0; i<npart; i++) {
          for (j=0; j<3; j++) {
            real_pos[i*3+j] = *(double*)((char*)PyArray_DATA(pos) + i*PyArray_STRIDES(pos)[0] + j*PyArray_STRIDES(pos)[1]);
          }
          real_value[i] = *(double*)((char*)PyArray_DATA(value) + i*PyArray_STRIDES(value)[0]);
        }

        cellsizex = bx / nx;
        cellsizey = by / ny;

        if(!cic && !tsc)
          {
            printf("Using nearest point assignment scheme\n");
            int cell;
            for(n=0; n<npart; n++)
              {
                i = floor((real_pos[n*3 + axis0] -cx + 0.5 * bx) /  cellsizex);
                j = floor((real_pos[n*3 + axis1] -cy + 0.5 * by) /  cellsizey);
                cell = j + i * ny;
                grid[ cell ] += real_value[ n ];
              }
          }
        else if(cic)
          {
            printf("Using CIC assignment scheme\n");
            int cell1, cell2, cell3, cell4;
            int k, l;
            for(n=0; n<npart; n++)
              {
                double pos_x = (real_pos[n*3 + axis0] -cx + 0.5 * bx) /  cellsizex - 0.5;
                double pos_y = (real_pos[n*3 + axis1] -cy + 0.5 * by) /  cellsizey - 0.5;

                i = floor(pos_x);
                j = floor(pos_y);

                pos_x -= i;
                pos_y -= j;

                k = j + 1;
                l = i + 1;
                /* periodic bcs */
                if(i < 0) 
                  i = nx - 1;
                if(j < 0) 
                  j = ny - 1;
                if(l >= nx)
                  l = 0;
                if(k >= ny)
                  k = 0;

                cell1 = j + i * ny;
                cell2 = k + i * ny;
                cell3 = j + l * ny;
                cell4 = k + l * ny;

                if(pos_x < 0.0 || pos_y < 0.0)
                  {
                    PyErr_SetString( PyExc_ValueError, "pos_x and pos_y must be positive" );
                    return 0;
                  }
                if(pos_x > 1.0 || pos_y > 1.0)
                  {
                    PyErr_SetString( PyExc_ValueError, "pos_x and pos_y must be <= 1" );
                    return 0;
                  }

/* non-periodic BCs
                if(i < 0) 
                  {
                    cell1 = j; 
                    cell2 = j + 1;
                    if(j < 0)
                      cell1 = cell2;
                    if(j >= ny - 1)
                      cell2 = cell1;
                  }
                if(j < 0) 
                  {
                    cell1 = i * ny; 
                    cell3 = (i + 1) * ny;
                    if(i < 0)
                      cell1 = cell3;
                    if(i >= nx - 1)
                      cell3 = cell1;
                  }

                if(i >= nx - 1) 
                  {
                    cell3 = j + (nx - 1) * ny; 
                    cell4 = (j + 1) + (nx - 1) * ny;
                    if(j < 0)
                      cell3 = cell4;
                    if(j >= ny - 1)
                      cell4 = cell3;
                  }
                if(j >= ny - 1) 
                  {
                    cell2 = (ny - 1) + i * ny; 
                    cell4 = (ny - 1) + (i + 1) * ny;
                    if(i < 0)
                      cell2 = cell4;
                    if(i >= nx - 1)
                      cell4 = cell2;
                  }
*/
                grid[ cell1 ] += real_value[ n ] * (1.0 - pos_y) * (1.0 - pos_x);
                grid[ cell2 ] += real_value[ n ] * pos_y * (1.0 - pos_x);
                grid[ cell3 ] += real_value[ n ] * (1.0 - pos_y) * pos_x;
                grid[ cell4 ] += real_value[ n ] * pos_y * pos_x;
              }
          }
        else if(tsc)
          {
            printf("Using TSC assignment scheme\n");
            PyErr_SetString( PyExc_ValueError, "TSC interpolation not implemented yet" );
            return 0;
          }

        free( real_pos );
        free( real_value );

        dict = PyDict_New();
        PyDict_SetStolenItem( dict, "grid", (PyObject*)pyGrid );
        
        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return dict;
}

PyObject* _calcDensProjection(PyObject *self, PyObject *args) {
        PyArrayObject *pos, *hsml, *mass, *pyGrid;
        int npart, nx, ny, cells;
        int dims[2];
        double bx, by, bz, cx, cy, cz;
        double *data_pos, *data_hsml, *data_mass;
        double *grid;
        int part;
        double px, py, pz, h, v, cpx, cpy, r2, h2;
        double p[3], weight;
        int x, y;
        int xmin, xmax, ymin, ymax, axis0, axis1;
        double cellsizex, cellsizey;
        clock_t start;
        
        start = clock();

        axis0 = 0;
        axis1 = 1;
        if (!PyArg_ParseTuple( args, "O!O!O!iidddddd|ii:calcDensProjection( pos, hsml, mass, nx, ny, boxx, boxy, boxz, centerx, centery, centerz, [axis0, axis1] )", &PyArray_Type, &pos, &PyArray_Type, &hsml, &PyArray_Type, &mass, &nx, &ny, &bx, &by, &bz, &cx, &cy, &cz, &axis0, &axis1 )) {
                return 0;
        }

        if (PyArray_NDIM(pos) != 2 || PyArray_DIMS(pos)[1] != 3 || PyArray_TYPE(pos) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "pos has to be of dimensions [n,3] and type double" );
                return 0;
        }

        if (PyArray_NDIM(hsml) != 1 || PyArray_TYPE(hsml) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "hsml has to be of dimension [n] and type double" );
                return 0;
        }

        if (PyArray_NDIM(mass) != 1 || PyArray_TYPE(mass) != NPY_DOUBLE) {
                PyErr_SetString( PyExc_ValueError, "mass has to be of dimension [n] and type double" );
                return 0;
        }

        npart = PyArray_DIMS(pos)[0];
        if (npart != PyArray_DIMS(hsml)[0] || npart != PyArray_DIMS(mass)[0]) {
                PyErr_SetString( PyExc_ValueError, "pos, hsml and mass have to have the same size in the first dimension" );
                return 0;
        }

        dims[0] = nx;
        dims[1] = ny;
        pyGrid = (PyArrayObject *)PyArray_FromDims( 2, dims, NPY_DOUBLE );
        grid = (double*)PyArray_DATA(pyGrid);
        cells = nx*ny;
        memset( grid, 0, cells*sizeof(double) );

        cellsizex = bx / nx;
        cellsizey = by / ny;

        data_pos = (double*)PyArray_DATA(pos);
        data_hsml = (double*)PyArray_DATA(hsml);
        data_mass = (double*)PyArray_DATA(mass);

        for (part=0; part<npart; part++) {
                p[0] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[1] = *data_pos;
                data_pos = (double*)((char*)data_pos + PyArray_STRIDES(pos)[1]);
                p[2] = *data_pos;
                data_pos = (double*)((char*)data_pos - 2*PyArray_STRIDES(pos)[1] + PyArray_STRIDES(pos)[0]);
                
                px = p[ axis0 ];
                py = p[ axis1 ];
                pz = p[ 3 - axis0 - axis1 ];
                
                h = *data_hsml;
                data_hsml = (double*)((char*)data_hsml + PyArray_STRIDES(hsml)[0]);
                h2 = h*h;

                v = *data_mass;
                data_mass = (double*)((char*)data_mass + PyArray_STRIDES(mass)[0]);

                xmin = max( floor( (px - h - cx + 0.5*bx) / cellsizex ), 0 );
                xmax = min( ceil( (px + h - cx + 0.5*bx) / cellsizex ), nx-1 );
                ymin = max( floor( (py - h - cy + 0.5*by) / cellsizey ), 0 );
                ymax = min( ceil( (py + h - cy + 0.5*by) / cellsizey ), ny-1 );
                
                if (xmin < nx && ymin < ny && xmax >= 0 && ymax >= 0 && pz > cz - 0.5*bz && pz < cz + 0.5*bz) {
                        weight = 0;
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        r2 = sqr(px-cpx-cx) + sqr(py-cpy-cy);
                                        if (r2 > h2) continue;
                                        weight += h * _getkernel( h, r2 );
                                }
                        }
                        
                        if (weight > 0)
                          {
                        for (x=xmin; x<=xmax; x++) {
                                cpx = -0.5*bx + bx*(x+0.5)/nx;
                                for (y=ymin; y<=ymax; y++) {
                                        cpy = -0.5*by + by*(y+0.5)/ny;
                                        r2 = sqr(px-cpx-cx) + sqr(py-cpy-cy);
                                        if (r2 > h2) continue;
                                        grid[x*ny + y] += h * _getkernel( h, r2 ) * v / weight;
                                }
                        }
          }
                }
        }

        printf( "Calculation took %gs\n", ((double)clock()-(double)start)/CLOCKS_PER_SEC );
        return PyArray_Return( pyGrid );
}

static PyMethodDef calcGridmethods[] = {
        { "calcGrid", (PyCFunction)_calcGrid, METH_VARARGS | METH_KEYWORDS, "" },
        { "calcSlice", _calcSlice, METH_VARARGS, "" },
        { "calcDensGrid", _calcDensGrid, METH_VARARGS, "" },
        { "calcDensSlice", _calcDensSlice, METH_VARARGS, "" },
        { "calcGridMassWeight", _calcGridMassWeight, METH_VARARGS, "" },
        { "calcCylinderAverage", _calcCylinderAverage, METH_VARARGS, "" },
        { "calcRadialProfile", _calcRadialProfile, METH_VARARGS, "" },
        { "calcAbundGrid", _calcAbundGrid, METH_VARARGS, "" },
        { "calcAbundGridSph", _calcAbundGridSph, METH_VARARGS, "" },
        { "calcAbundSphere", _calcAbundSphere, METH_VARARGS, "" },
        { "gatherAbundGrid", (PyCFunction)_gatherAbundGrid, METH_VARARGS | METH_KEYWORDS, "" },
        { "calcASlice", (PyCFunction)_calcASlice, METH_VARARGS | METH_KEYWORDS, "" },
        { "calcASlice2", (PyCFunction)_calcASlice2, METH_VARARGS | METH_KEYWORDS, "" },
        { "calcDensProjection", _calcDensProjection, METH_VARARGS, "" },
        { NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "calcGrid", /* m_name */
  NULL,      /* m_doc */
  -1,                  /* m_size */
  calcGridmethods,     /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_calcGrid(void)
{
        import_array();
        return PyModule_Create( &moduledef );
}
#else
PyMODINIT_FUNC initcalcGrid(void)
{
        Py_InitModule( "calcGrid", calcGridmethods );
        import_array();
}
#endif
