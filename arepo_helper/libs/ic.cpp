#include <Python.h>
#include <arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "helm_eos.h"
#include "pyhelm_eos.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "ic.h"
#include "const.h"
#include "utils.h"
#include "spline.h"

int create_wd_integrator(double r, const double *y, double *ydot, void *params)
{
    double temp = ((struct paramsWD*)params)->temp;
    double *xnuc = ((struct paramsWD*)params)->xnuc;
    t_helm_eos_table *eos = ((struct paramsWD*)params)->eos;
    double rho = ((struct paramsWD*)params)->rho;

    struct eos_result res{};
    eos_calc_ptgiven(eos, y[0], xnuc, temp, &rho, &res);
    r = std::max(r, 1e2);
    double r2 = r * r;

    ydot[0] = - G * y[1] * rho / r2;
    ydot[1] = 4. * M_PI  * r2 * rho;

    ((struct paramsWD*)params)->rho = rho;

    return GSL_SUCCESS;
}

class WdecResults {

public:
    int n_points{};
    std::string output_filepath;
    std::string corsico_filepath;
    std::vector<double> radii;
    std::vector<double> mr;
    std::vector<double> temp;
    std::vector<double> density;
    std::vector<double> pres;
    std::vector<double> q;
    std::vector<double> xo;
    std::vector<double> xc;
    std::vector<double> xhe;
    std::vector<double> xh;


    explicit WdecResults(const std::string& directory) {

        output_filepath     = directory + "output.dat";
        corsico_filepath    = directory + "corsico.dat";

        read_output_dat();
        read_corsico_dat();

        n_points -= 2;
        truncate_vector(&radii);
        truncate_vector(&mr);
        truncate_vector(&temp);
        truncate_vector(&density);
        truncate_vector(&pres);

        // Fill in the centre of the star
        int num_infill      = 0;
        radii = extend_vector(&radii, 100, radii[0] - 1000, num_infill);
        mr = extend_vector(&mr, 100, mr[0] - 1000, num_infill);
        temp = extend_vector(&temp, temp[0], temp[0], num_infill);
        density = extend_vector(&density, density[0], density[0], num_infill);
        pres = extend_vector(&pres, pres[0], pres[0], num_infill);
        xo = extend_vector(&xo, xo[0], xo[0], num_infill);
        xh = extend_vector(&xh, xh[0], xh[0], num_infill);
        xhe = extend_vector(&xhe, xhe[0], xhe[0], num_infill);
        xc = extend_vector(&xc, xc[0], xc[0], num_infill);
        n_points += num_infill;

        xc.pop_back();
        xc.pop_back();
        xo.pop_back();
        xo.pop_back();
    };

    void read_output_dat() {

        std::string line;
        std::ifstream file;
        std::stringstream ss;
        file.open(output_filepath, std::fstream::binary | std::fstream::out);
        if (!file)
            std::cout << "File Not Found." << std::endl;
        else {
            file.seekg(0); // To make sure that the data is read from the starting position of the file.

            // skip 2 lines
            std::getline(file, line);
            std::getline(file, line);

            // "npoints"
            std::getline(file, line);
            std::stringstream ss(line);
            ss >> n_points;

            // skip 2 lines
            std::getline(file, line);
            std::getline(file, line);

            for (int i = 0; i < n_points; i++) {
                std::getline(file, line);
                std::stringstream ss(line);
                double ni, ri, mri, lri, tempi, rhoi, presi;
                ss >> ni >> ri >> mri >> lri >> tempi >> rhoi >> presi;

                radii.push_back(ri);
                mr.push_back(mri);
                temp.push_back(tempi);
                density.push_back(rhoi);
                pres.push_back(presi);

            }

            file.close();
        }
    }

    void read_corsico_dat() {
        std::string line;
        std::ifstream file;
        file.open(corsico_filepath, std::fstream::binary | std::fstream::out);
        if (!file)
            std::cout << "File Not Found." << std::endl;
        else {
            file.seekg(0); // To make sure that the data is read from the starting position of the file.

            // skip 1 lines
            std::getline(file, line);


            for (int i = 0; i < n_points; i++) {
                std::getline(file, line);
                std::stringstream ss(line);
                double qi, xoi, xci, xhei, xhi;
                ss >> qi >> xoi >> xci >> xhei >> xhi;

                q.push_back(qi);
                xo.push_back(xoi);
                xc.push_back(xci);
                xhe.push_back(xhei);
                xh.push_back(xhi);

            }

            file.close();
        }
    }

    static void truncate_vector(std::vector<double> *a) {
        a->pop_back();
        a->erase(a->begin());
    }

    template<typename T> static std::vector<double> linspace(T start_in, T end_in, int num_in) {

        std::vector<double> linspaced;

        auto start = static_cast<double>(start_in);
        auto end = static_cast<double>(end_in);
        auto num = static_cast<double>(num_in);

        if (num == 0) { return linspaced; }
        if (num == 1)
        {
            linspaced.push_back(start);
            return linspaced;
        }

        double delta = (end - start) / (num - 1);

        for(int i=0; i < num-1; ++i)
        {
            linspaced.push_back(start + delta * i);
        }
        linspaced.push_back(end); // I want to ensure that start and end
        // are exactly the same as the input
        return linspaced;
    }

    static std::vector<double> extend_vector(std::vector<double> *vec, double start, double end, int num_infill) {
        std::vector<double> extended = linspace(start, end, num_infill);
        extended.insert(extended.end(), vec->begin(), vec->end());

        return extended;
    }
};

double mtot_from_rho_c_implementation(double rho_c,
                                      double temp_c,
                                      t_helm_eos_table *eos,
                                      PyObject *xnuc_py,
                                      double offset) {
    auto wd_guess = create_wd_implementation(eos, rho_c, temp_c, xnuc_py, 1e-6);
    auto mr_guess = (PyArrayObject *) PyDict_GetItemString(wd_guess, f[N::MR]);
    auto mr_data = (double *) PyArray_DATA(mr_guess);
    auto npart = PyArray_DIMS(mr_guess)[0];

    return mr_data[npart - 1] / msol - offset;
}

double rho_c_from_mtot_implementation(double mtot,
                                      double temp_c,
                                      t_helm_eos_table *eos,
                                      PyObject *xnuc_py) {

    double tolerance    = 0.01;
    double rho_guess    = 1e4;
    double mass_guess   = mtot_from_rho_c_implementation(rho_guess, temp_c, eos, xnuc_py, mtot);

    while (mass_guess < 0) {
        rho_guess   *= 10.0;
        mass_guess  = mtot_from_rho_c_implementation(rho_guess, temp_c, eos, xnuc_py, mtot);
    }

    auto a = 1.0/10.0 * rho_guess;
    auto b = rho_guess;

    std::cout << "Finding mass in range [a, b] = " << a << ", " << b << std::endl;

    double lower, mid;
    while ((b - a) > tolerance) {

        // Find middle point
        rho_guess = (a + b) / 2;
        mid = mtot_from_rho_c_implementation(rho_guess, temp_c, eos, xnuc_py, mtot);
        lower = mtot_from_rho_c_implementation(a, temp_c, eos, xnuc_py, mtot);

        // Check if middle point is root
        if (mid == 0.0)
            break;

        // Else, shrink the search space
        else if (mid * lower < 0)
            b = rho_guess;
        else
            a = rho_guess;
    }

    std::cout << "Found rho_c = " << rho_guess << " . Mass [MSol] = " <<
                                                                      mtot_from_rho_c_implementation(rho_guess, temp_c,
                                                                                                     eos, xnuc_py) << std::endl;

    return rho_guess;
}

PyObject *create_wd_implementation(t_helm_eos_table *eos,
                                   double rho_c,
                                   double temp_c,
                                   PyObject *xnuc_py,
                                   double tolerance) {


    int max_profile_length  = 0x7FFFFFFF;
    double rho_cutoff       = 1e-5;
    double r_cutoff         = 1e10;
    double r_c              = 1e2;
    double r_c_deriv        = 1e3;
    int count               = 0;

    int nspecies = PyList_Size(xnuc_py);
    auto xnuc = (double *) malloc(nspecies * sizeof(double));
    for (int i = 0; i < nspecies; i++) {
        xnuc[i] = PyFloat_AsDouble(PyList_GetItem(xnuc_py, i));
    }


    // Initialise the profile primitive variable containers
    auto r = (double *) malloc(max_profile_length * sizeof(double));
    auto p = (double *) malloc(max_profile_length * sizeof(double));
    auto u = (double *) malloc(max_profile_length * sizeof(double));
    auto rho = (double *) malloc(max_profile_length * sizeof(double));
    auto dm = (double *) malloc(max_profile_length * sizeof(double));
    auto mr = (double *) malloc(max_profile_length * sizeof(double));
    auto csnd = (double *) malloc(max_profile_length * sizeof(double));
    auto xnuc_radial = (double *) malloc(nspecies * max_profile_length * sizeof(double));

    // Find the state in the centre
    struct eos_result res{};
    eos_calc_tgiven(eos, rho_c, xnuc, temp_c, &res);
    double initial_mass = rho_c * 4.0 / 3.0 * M_PI * r_c * r_c * r_c;
    dm[0]               = initial_mass;
    mr[0]               = initial_mass;
    rho[0]              = rho_c;
    r[0]                = r_c;
    p[0]                = res.p.v;
    u[0]                = res.e.v;
    csnd[0]             = res.sound;
    for (int i = 0; i < nspecies; i++) {
        xnuc_radial[0 + i] = xnuc[i];
    }

    count++;

    // Initialise the ODE solver
    struct paramsWD params {.temp = temp_c, .xnuc = xnuc, .eos = eos, .rho = rho_c};
    const gsl_odeiv_step_type * T   = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * s              = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c           = gsl_odeiv_control_y_new (0.0, tolerance);
    gsl_odeiv_evolve * ev           = gsl_odeiv_evolve_alloc (2);
    gsl_odeiv_system sys            = {create_wd_integrator, nullptr, 2, &params};

    // Pressure and total mass within radius r
    double y[2]     = {p[0], 0.0};
    double rad      = r_c;
    double drad     = r_c_deriv;
    bool success    = true;

    while (rho[count - 1] > rho_cutoff && rad < r_cutoff && count < max_profile_length) {

        int evolution_status = gsl_odeiv_evolve_apply(ev, c, s, &sys, &rad, 1e10, &drad, y);

        // Quit the loop early if evolution failed
        if (evolution_status != GSL_SUCCESS || y[0] < 0.0) {
            success = false;
            break;
        }

        // Compute the new values of the primitive quantities
        rho[count]  = rho[count-1];
        eos_calc_ptgiven(eos, y[0], xnuc, temp_c, &rho[count], &res);
        r[count]    = rad;
        p[count]    = y[0];
        u[count]    = res.e.v;
        csnd[count] = res.sound;
        mr[count]   = y[1];
        dm[count]   = y[1] - mr[count - 1];

        for (int i = 0; i < nspecies; i++) {
            xnuc_radial[count * nspecies + i] = xnuc[i];
        }

        count++;
    }

    // Free all unused space
    r = (double *) realloc(r, count * sizeof(double));
    p = (double *) realloc(p, count * sizeof(double));
    u = (double *) realloc(u, count * sizeof(double));
    rho = (double *) realloc(rho, count * sizeof(double));
    dm = (double *) realloc(dm, count * sizeof(double));
    mr = (double *) realloc(mr, count * sizeof(double));
    csnd = (double *) realloc(csnd, count * sizeof(double));
    xnuc_radial = (double *) realloc(xnuc_radial, nspecies * count * sizeof(double));

    PyObject* dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::DENSITY], (PyObject *) create_numpy_array(rho, count));
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject *) create_numpy_array(u, count));
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject *) create_numpy_array(p, count));
    PyDict_SetStolenItem(dict, f[N::RADIUS], (PyObject *) create_numpy_array(r, count));
    PyDict_SetStolenItem(dict, f[N::DM], (PyObject *) create_numpy_array(dm, count));
    PyDict_SetStolenItem(dict, f[N::MR], (PyObject *) create_numpy_array(mr, count));
    PyDict_SetStolenItem(dict, f[N::SOUNDSPEED], (PyObject *) create_numpy_array(csnd, count));
    PyDict_SetStolenItem(dict, f[N::NUCLEARCOMPOSITION], (PyObject *) create_numpy_array(xnuc_radial, count, nspecies));

    free(r);
    free(p);
    free(u);
    free(rho);
    free(dm);
    free(mr);
    free(csnd);
    free(xnuc_radial);
    free(xnuc);

    gsl_odeiv_evolve_free(ev);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);

    if (success) {
        return dict;
    } else {
        return nullptr;
    }

}

PyObject *create_polytrope_implementation(t_helm_eos_table *eos,
                                          double n,
                                          double rho_c,
                                          PyObject *xnuc_py,
                                          double pres_c,
                                          double temp_c,
                                          double dr) {

    int max_profile_length  = 0x7FFFFFFF;
    double density_cutoff   = 1.0e-3;
    auto xnuc               = (double *) xnuc_py;

    // Allocate memory for profiles
    auto r      = (double *) malloc(max_profile_length * sizeof(double));
    auto rho    = (double *) malloc(max_profile_length * sizeof(double));
    auto p      = (double *) malloc(max_profile_length * sizeof(double));
    auto u      = (double *) malloc(max_profile_length * sizeof(double));
    auto dm     = (double *) malloc(max_profile_length * sizeof(double));
    auto mr     = (double *) malloc(max_profile_length * sizeof(double));

    // Find the state in the center
    struct eos_result res{};
    double temp;

    if (pres_c == 0.0) {
        temp = temp_c;
        eos_calc_tgiven(eos, rho_c, xnuc, temp, &res);
        p[0] = res.p.v;
    } else {
        temp = -1.0;
        eos_calc_pgiven(eos, rho_c, xnuc, pres_c, &temp, &res);
        p[0] = pres_c;
    }

    rho[0]  = rho_c;
    u[0]    = res.e.v;
    dm[0]   = 0.0;
    r[0]    = 0.0;

    // Find the power law constant K and display info
    double gamma            = (n + 1.0) / n;
    double one_over_gamma   = 1.0 / gamma;
    double K                = p[0] / pow(rho[0], gamma);
    printf("Central density: %g\n", rho[0]);
    printf("Central temperature: %g\n", temp);
    printf("Polytropic index: %g\n", n);
    printf("Polytropic constant: %g\n", K);

    // Fill in further profile points using dpdr and eos_calc_pgiven
    int i = 1;
    double mass = 0;
    double rad = 0;
    double dpdr, density_av, shell_volume, r_squared;

    while (i < max_profile_length && rho[i - 1] > density_cutoff) {

        rad         += dr;
        r[i]        = rad;
        r_squared   = std::max(r[i - 1] * r[i - 1], 0.1);

        dpdr        = -G * mass * rho[i - 1] / r_squared;
        p[i]        = p[i-1] + dr * dpdr;
        rho[i]      = pow(p[i] / K, one_over_gamma);


        eos_calc_pgiven(eos, rho[i], xnuc, p[i], &temp, &res);
        u[i] = res.e.v;

        density_av      = 0.5 * (rho[i-1] + rho[i]);
        shell_volume    = 4.0 * M_PI * pow((rad - 0.5 * dr), 2.0) * dr;
        dm[i] = density_av * shell_volume;
        mass += dm[i];
        mr[i] = mass;

        i++;
    }

    // Free unused memory
    r      = (double *) realloc(r, i * sizeof(double));
    rho    = (double *) realloc(rho, i * sizeof(double));
    p      = (double *) realloc(p, i * sizeof(double));
    u      = (double *) realloc(u, i * sizeof(double));
    dm     = (double *) realloc(dm, i * sizeof(double));
    mr     = (double *) realloc(mr, i * sizeof(double));

    // Write profiles to Python dict and free all other memory
    auto dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::RADIUS], (PyObject*) create_numpy_array(r, i));
    PyDict_SetStolenItem(dict, f[N::DENSITY], (PyObject*) create_numpy_array(rho, i));
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject*) create_numpy_array(u, i));
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject*) create_numpy_array(p, i));
    PyDict_SetStolenItem(dict, f[N::DM], (PyObject*) create_numpy_array(dm, i));
    PyDict_SetStolenItem(dict, f[N::MR], (PyObject*) create_numpy_array(mr, i));

    free(r);
    free(rho);
    free(dm);
    free(mr);
    free(p);
    free(u);

    return dict;
}

PyObject *create_wd_wdec_implementation(const char *wdec_dir, int nspecies, double gamma) {

    auto wd_results = WdecResults(wdec_dir);

    int npoints = wd_results.n_points;
    auto r     = convert_to_double_star(wd_results.radii, 1);
    auto mr_old    = convert_to_double_star(wd_results.mr, 1);
    auto temp  = convert_to_double_star(wd_results.temp, 1);
    auto q      = convert_to_double_star(wd_results.q, 1);
    auto rho   = convert_to_double_star(wd_results.density, 1);
    auto p     = convert_to_double_star(wd_results.pres, 1);
    auto xo    = convert_to_double_star(wd_results.xo, 1);
    auto xc    = convert_to_double_star(wd_results.xc, 1);

    auto xnuc_radial = (double *) malloc(nspecies * npoints * sizeof(double));
    double xCHere, xOHere;
    std::vector<double> u_vec, dm_vec, mr_vec;
    double delta_mr, central_mass;
    u_vec.reserve(npoints);
    dm_vec.reserve(npoints);
    mr_vec.reserve(npoints);

    for (int i = 0; i < npoints; i++) {
        // Fill in energy
        u_vec.push_back(p[i] / (rho[i] * (gamma - 1)));

        // Fill in mass difference between shells
        if (i == 0) {
            central_mass = 4.0 / 3.0 * M_PI * rho[0] * r[0] * r[0] * r[0];
            dm_vec.push_back(central_mass);
            mr_vec.push_back(central_mass);
        } else {
            delta_mr = mr_diff(q[i - 1], q[i]);
            dm_vec.push_back(delta_mr);
            mr_vec.push_back(mr_vec[i - 1] + delta_mr);
        }

        double mr_diff = mr_vec[i] - mr_old[i];
//        std::cout << mr_diff << std::endl;

        // Fill in nuclear compositions
        for (int j = 0; j < nspecies; j++) {
            xnuc_radial[i * nspecies + j] = 0.0;
        }

        // Specifically fill in C, O, He and H
        xOHere = xo[i];
        xCHere = xc[i];

        if (xCHere + xOHere >= 1) {xOHere = 1 - xCHere;}

        xnuc_radial[i * nspecies + 1] = xCHere; // C
        xnuc_radial[i * nspecies + 2] = xOHere; // O
        xnuc_radial[i * nspecies + 0] = 1.0 - xCHere - xOHere; // He
    }

    auto u = convert_to_double_star(u_vec, 1);

    // Kludge
    double mr_true = mr_old[npoints - 1];
    double mr_fake = mr_vec[npoints - 1];
    for (int i = 0; i < npoints; i++) {
        mr_vec[i] *= mr_true / mr_fake;
        dm_vec[i] *= mr_true / mr_fake;
    }

    auto dm = convert_to_double_star(dm_vec, 1);
    auto mr = convert_to_double_star(mr_vec, 1);

    double mtot = mr[npoints - 1];
    printf("Total mass: %g solar masses\n", mtot / msol);


    auto dict = PyDict_New();
    PyDict_SetStolenItem(dict, f[N::RADIUS], (PyObject*) create_numpy_array(r, npoints));
    PyDict_SetStolenItem(dict, f[N::DENSITY], (PyObject*) create_numpy_array(rho, npoints));
    PyDict_SetStolenItem(dict, f[N::TEMPERATURE], (PyObject*) create_numpy_array(temp, npoints));
    PyDict_SetStolenItem(dict, f[N::INTERNALENERGY], (PyObject*) create_numpy_array(u, npoints));
    PyDict_SetStolenItem(dict, f[N::PRESSURE], (PyObject*) create_numpy_array(p, npoints));
    PyDict_SetStolenItem(dict, f[N::DM], (PyObject*) create_numpy_array(dm, npoints));
    PyDict_SetStolenItem(dict, f[N::MR], (PyObject*) create_numpy_array(mr, npoints));
    PyDict_SetStolenItem(dict, f[N::NUCLEARCOMPOSITION], (PyObject*) create_numpy_array(xnuc_radial, npoints, nspecies));

    return dict;

}

PyObject *create_wd(PyObject *self, PyObject *args, PyObject *kwargs) {
    t_helm_eos_table *eos;
    double rho_c;
    double temp_c           = 5e5;
    double tolerance        = 1e-6;
    double xnuc_tolerance   = 1e-14;
    double xtot             = 0.0;
    PyObject *xnuc_py;

    const char *kwlist[] = {"eos", "rho_c", "temp_c", "xnuc_py", "tolerance", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O&d|dOd:create_wd(eos, rho_c, [temp_c, xnuc_py, tolerance])",
                                     keywords,
                                     &pyConvertHelmEos, &eos, &rho_c, &temp_c, &xnuc_py, &tolerance)) {
        return nullptr;
    }

    if (PyList_Size(xnuc_py) != eos->nspecies) {
        PyErr_SetString(PyExc_ValueError, "Nspecies mismatch between EOS and xnuc_py. \n");
        return nullptr;
    }

    for (int i = 0; i < eos->nspecies; i++) {
        xtot += PyFloat_AsDouble(PyList_GetItem(xnuc_py, i));
    }

    if(fabs(xtot - 1.0) > xnuc_tolerance)
    {
        PyErr_SetString(PyExc_ValueError, "Inconsistent Abundances. \n");
        return nullptr;
    }

    auto dict = create_wd_implementation(eos, rho_c, temp_c, xnuc_py, tolerance);

    if (dict == nullptr) {
    PyErr_SetString(PyExc_ValueError, "ODE Solver failed.\n");
    return nullptr;
    }


    return dict;
}

PyObject *create_polytrope(PyObject *self, PyObject *args, PyObject *kwargs) {
    t_helm_eos_table *eos;
    double rho_c;
    double n;
    double pres_c = 0.0;
    double temp_c = 0.0;
    double dr = 1e5;
    PyObject *xnuc_py;

    const char *kwlist[] = {"eos", "n", "rho_c", "xnuc_py", "pres_c", "temp_c", "dr", nullptr};
    auto keywords = (char **) kwlist;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                     "O&ddO|ddd:create_polytrope(eos, n, rho_c, xnuc_py, [pres_c, temp_c, dr])",
                                     keywords,
                                     &pyConvertHelmEos, &eos,
                                     &n, &rho_c,
                                     &xnuc_py,
                                     &pres_c, &temp_c, &dr)) {
        return nullptr;
    }

    if (PyList_Size(xnuc_py) != eos->nspecies) {
        PyErr_SetString(PyExc_ValueError, "Nspecies mismatch between EOS and xnuc_py. \n");
        return nullptr;
    }

    if (pres_c == 0.0 && temp_c == 0.0) {
        PyErr_SetString(PyExc_ValueError, "pres_c or temp_c have to be set.\n");
        return nullptr;
    }

    auto dict = create_polytrope_implementation(eos, n, rho_c, xnuc_py, pres_c, temp_c, dr);

    return dict;
}

PyObject *create_wd_wdec(PyObject *self, PyObject *args, PyObject *kwargs) {
    char *wdec_dir;
    int nspecies;
    double gamma = 1.4;

    const char *kwlist[] = {"wdec_dir", "nspecies", "gamma", nullptr};
    auto keywords = (char **) kwlist;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "si|d:create_wd_wdec(wdec_dir, nspecies, [gamma])",
                          keywords, &wdec_dir, &nspecies, &gamma)) {
        return nullptr;
    }
    auto dict = create_wd_wdec_implementation(wdec_dir, nspecies, gamma);
    return dict;
}

PyObject *mtot_from_rho_c(PyObject *self, PyObject *args) {

    double rho_c;
    double temp_c;
    t_helm_eos_table *eos;
    PyObject *xnuc_py;

    if (!PyArg_ParseTuple(args, "ddO&O:mtot_from_rho_c()",
                          &rho_c, &temp_c,
                          &pyConvertHelmEos, &eos,
                          &xnuc_py)) {
        PyErr_SetString(PyExc_ValueError, "Parse error!");
        return nullptr;
    }

    auto mtot = mtot_from_rho_c_implementation(rho_c, temp_c, eos, xnuc_py);

    return PyFloat_FromDouble(mtot);
}

PyObject *rho_c_from_mtot(PyObject *self, PyObject *args) {
    double mtot;
    double temp_c;
    t_helm_eos_table *eos;
    PyObject *xnuc_py;

    if (!PyArg_ParseTuple(args, "ddO&O:rho_c_from_mtot()",
                          &mtot, &temp_c,
                          &pyConvertHelmEos, &eos,
                          &xnuc_py)) {
        PyErr_SetString(PyExc_ValueError, "Parse error!");
        return nullptr;
    }

    auto rho_c = rho_c_from_mtot_implementation(mtot, temp_c, eos, xnuc_py);

    return PyFloat_FromDouble(rho_c);
}

PyMODINIT_FUNC PyInit_ic(void)
{
    import_array()
    return PyModule_Create(&moduledef_ic);
}
