#ifndef AREPO_HELPER_LIBS_MAKE_WDEC_H
#define AREPO_HELPER_LIBS_MAKE_WDEC_H

#include <Python.h>
#include <hdf5.h>
#include <arrayobject.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

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
        int num_infill      = 200;
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

    template<typename T>
    static std::vector<double> linspace(T start_in, T end_in, int num_in) {

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

void print_vector(const std::vector<double>& vec);
//inline PyArrayObject *createPyArray(double *data, int dim1, int dim2);
std::vector<double> extend_vector(std::vector<double> *vec, double start, double end, int num_infill);
double *convert_to_double_star(std::vector<double> vector, int n);
void write_the_hdf5_stuff(hid_t hdf5_file,
                          hid_t hdf5_group,
                          int num_particles,
                          const std::vector<const char *>& field_names,
                          const std::vector<int>& num_columns_v,
                          const std::vector<double *>& all_data);
void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, double value);
void write_scalar_to_file(hid_t headergrp, const char *name, hid_t datatype, int value);
void write_vector_to_file(hid_t headergrp, int ntypes, const char *field_name, hid_t datatype, int data[6]);

static void write_float_buffer_vector(hid_t file_id,
                          unsigned int nof_particles,
                          double *buffer,
                          const char *dataset_name,
                          long hdf_datatype,
                          int num_columns);

static void write_float_buffer(hid_t file_id,
                               unsigned int nof_particles,
                               double *buffer,
                               const char *dataset_name,
                               long hdf_datatype);


PyObject *make_wdec_new(const char *wdec_dir,
                   double boxsize,
                   int nspecies,
                   bool makebox,
                   bool randomizeshells,
                   bool randomizeradii,
                   double pmass);

#endif //AREPO_HELPER_LIBS_MAKE_WDEC_H
