//
//  NonGaussianParameter.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "NonGaussianParameter.hpp"

#ifdef OMP
#include "omp.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <iomanip>
#include <cstring>

#include "Trajectory.hpp"

using namespace std;

NonGaussianParameter::NonGaussianParameter()
{
    input_file_name_ = "alpha2_t.in";
    output_file_name_ = "alpha2_t.txt";
}


NonGaussianParameter::~NonGaussianParameter()
{
}


void NonGaussianParameter::set_output_file_name()
{
    if (output_file_name_ == "r2_t.txt") {
        output_file_name_ = "alpha2_t.txt";
    }
}


void NonGaussianParameter::compute_alpha2_t()
{
    if (is_wrapped_) {
        unwrap_coordinates();
    }
    
    r2_t_.resize(number_of_time_points_, 0.0);
    alpha2_t_.resize(number_of_time_points_, 0.0);
    
    // Form Array of time index values for a given type of timescale computation
    compute_time_array();
    
    // select the indexes of atom_type_ in atom_group_
    vector< unsigned int > atom_type_indexes;
    select_atoms(atom_type_indexes, atom_type_, atom_group_);
    
    double const normalization_factor = 1.0/(atom_type_indexes.size() * number_of_frames_to_average_);
    
    cout << setiosflags(ios::fixed);
    cout << setprecision(4);
    
    int status = 0;
    cout << "Computing ..." << endl;
    
    // Perform time averaging of Non Gaussian Parameter
#pragma omp parallel for
    for (size_t time_point = 1; time_point <  number_of_time_points_; ++time_point) {
        double total_squared_displacement = 0.0;
        double total_bisquared_displacement = 0.0;
        for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
            size_t current_frame = initial_frame + time_array_indexes_[time_point];
            for (size_t i_atom = 0; i_atom < atom_type_indexes.size(); ++i_atom) {
                size_t atom_index = atom_type_indexes[i_atom];
                double r_squared_dimension_sum = 0.0;
                for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                    double delta_x = trajectory_[current_frame][atom_index][i_dimension] - trajectory_[initial_frame][atom_index][i_dimension];
                    r_squared_dimension_sum += delta_x * delta_x;
                }
                total_squared_displacement += r_squared_dimension_sum;
                total_bisquared_displacement += r_squared_dimension_sum * r_squared_dimension_sum;
            }
        }
        r2_t_[time_point] = total_squared_displacement * normalization_factor;
        alpha2_t_[time_point] = 3.0 * total_bisquared_displacement/ (5.0 * total_squared_displacement * total_squared_displacement * normalization_factor) - 1.0;
        
        if (is_run_mode_verbose_) {
#pragma omp critical
            {
                ++status;
                cout << "\rcurrent progress of calculating non-gaussian parameter is: ";
                cout << status * 100.0/number_of_time_points_;
                cout << " \%";
                cout << flush;
            }
        }
    }
    
    cout << endl;
}


void NonGaussianParameter::write_alpha2_t()
{
    if (trajectory_delta_time_ == 0.0) {
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
    ofstream output_alpha2_t_file(output_file_name_);
    
    if (!output_alpha2_t_file) {
        cerr << "ERROR: Output file for non gaussian parameter: ";
        cerr << "\033[1;25m";
        cerr << output_file_name_;
        cerr << ", could not be opened.";
        cerr << endl;
        exit(1);
    }
    
    output_alpha2_t_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    output_alpha2_t_file << "#Non Gaussian Parameter and Mean Squared Displacement for ";
    output_alpha2_t_file << atom_type_;
    output_alpha2_t_file << " atoms of group ";
    output_alpha2_t_file << atom_group_;
    output_alpha2_t_file << "\n";
    output_alpha2_t_file << "#using ";
    output_alpha2_t_file << time_scale_type_;
    output_alpha2_t_file << "scale\n";
    output_alpha2_t_file << "#time               MSD                 NGP \n";
    
    for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
        output_alpha2_t_file << time_array_indexes_[time_point] * trajectory_delta_time_;
        output_alpha2_t_file << "        ";
        output_alpha2_t_file << r2_t_[time_point];
        output_alpha2_t_file << "        ";
        output_alpha2_t_file << alpha2_t_[time_point];
        output_alpha2_t_file << "\n";
    }

    output_alpha2_t_file.close();
}