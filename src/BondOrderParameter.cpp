//
//  BondOrderParameter.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "BondOrderParameter.hpp"

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
#include <random>
#include <sstream>
#include <cstring>

#ifdef GSL
#include <gsl/gsl_sf_legendre.h>
#endif

#ifdef BOOST
#include <boost/math/special_functions/spherical_harmonic.hpp>
#endif

#include "Trajectory.hpp"

using namespace std;

BondOrderParameter::BondOrderParameter() :
    is_averaged_(false),
    add_bar_(false),
    input_file_name_("BOP.in"),
    output_file_name_("BOP.txt"),
    time_scale_type_("linear"),
    number_of_time_points_(0),
    bond_parameter_order_(0),
    frame_interval_(1.0),
    max_cutoff_length_(0.0),
    atom_types_({}),
    atom_group_("system")
{
}


BondOrderParameter::~BondOrderParameter()
{
}


void BondOrderParameter::read_command_inputs(int argc, char * argv[])
{
    for (int input = 1; input < argc; ++input) {
        if (strcmp(argv[input], "-i") == 0) {
            input_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-o") == 0) {
            output_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-t") == 0) {
            trajectory_file_name_ = argv[++input];
            continue;
        }
        if (strcmp(argv[input], "-v") == 0) {
            is_run_mode_verbose_ = 1;
            continue;
        }
        
        cerr << "\nERROR: Unrecognized flag '" << argv[input] << "' from command inputs.\n";
        exit(1);
    }
}


void BondOrderParameter::read_input_file()
{
    ifstream input_file(input_file_name_);
    
    if (!input_file) {
        cerr << "ERROR: Input file location, ";
        cerr << "\033[1;25m";
        cerr << input_file_name_;
        cerr << "\033[0m";
        cerr << ", does not exist, please check input";
        cerr << endl;
        exit(1);
    }
    
    string input_word;
    
    while (input_file >> input_word) {
        
        //check for comment
        if (input_word[0] == '#') {
            getline(input_file, input_word);
            continue;
        }
        
        //check for memeber bools
        if (input_word == "is_averaged") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                is_averaged_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                is_averaged_ = false;
            }
            else {
                is_averaged_ = stoi(input_word);
            }
            continue;
        }
        if (input_word == "add_bar") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                add_bar_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                add_bar_ = false;
            }
            else {
                add_bar_ = stoi(input_word);
            }
            continue;
        }
        if (input_word == "is_run_mode_verbose") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                is_run_mode_verbose_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                is_run_mode_verbose_ = false;
            }
            else {
                is_run_mode_verbose_ = stoi(input_word);
            }
            continue;
        }
        if (input_word == "is_wrapped") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            if (input_word == "true" || input_word == "yes") {
                is_wrapped_ = true;
            }
            else if(input_word == "false" || input_word == "no") {
                is_wrapped_ = false;
            }
            else {
                is_wrapped_ = stoi(input_word);
            }
            continue;
        }
        
        //check if equal to member strings
        if (input_word == "output_file_name") {
            if (output_file_name_ != "BOP.txt") {
                cerr << "ERROR: Please do not set output file by command line and input file,\n";
                cerr << "     : we are unsure on which to prioritize.";
                cerr << endl;
                exit(1);
            }
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            output_file_name_ = input_word;
            continue;
        }
        if (input_word == "trajectory_file_name") {
            if (trajectory_file_name_ != "") {
                cerr << "ERROR: Please do not set trajectory file by command line and input file,\n";
                cerr << "     : we are unsure on which to prioritize.";
                cerr << endl;
                exit(1);
            }
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_file_name_ = input_word;
            continue;
        }
        if (input_word == "trajectory_file_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_file_type_ = input_word;
            continue;
        }
#ifdef GROMACS
        if (input_word == "gro_file_name") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            gro_file_name_ = input_word;
            continue;
        }
#else
        if (input_word == "gro_file_name") {
            cerr << "WARNING: gro files cannot be used in non gromacs";
            cerr << "compatible version of LiquidLib\n";
            cerr << endl;
        }
#endif
        
        if (input_word == "atom_group") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            atom_group_ = input_word;
            continue;
        }
        // Read atom types and scattering lengths provided by user
        if (input_word == "atom_types") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                atom_types_.push_back(input_word);
            }
            continue;
        }

        if (input_word == "time_scale_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            time_scale_type_ = input_word;
            continue;
        }
        if (input_word == "trajectory_data_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_data_type_ = input_word;
            continue;
        }
        
        //check if equal to member ints
        if (input_word == "start_frame") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            start_frame_ = stoi(input_word);
            continue;
        }
        if (input_word == "end_frame") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            end_frame_ = stoi(input_word);
            continue;
        }
        if (input_word == "dimension") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            dimension_ = stoi(input_word);
            continue;
        }
        if (input_word ==  "number_of_time_points") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            number_of_time_points_ = stoi(input_word);
            continue;
        }
        if (input_word ==  "bond_parameter_order") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            bond_parameter_order_ = stoi(input_word);
            continue;
        }
        
        //check if equal to member doubles
        if (input_word == "frame_interval") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            frame_interval_ = stod(input_word);
            continue;
        }
        if (input_word == "max_cutoff_length") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            max_cutoff_length_ = stod(input_word);
            continue;
        }
        if (input_word == "trajectory_delta_time") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_delta_time_ = stod(input_word);
            continue;
        }
        if (input_word == "output_precision") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            output_precision_ = stod(input_word);
            continue;
        }
        if (input_word == "boxlength") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                average_box_length_.push_back(stod(input_word));
            }
            continue;
        }
        
        //check for everything else
        cerr << "WARNING: no matching input type for: ";
        cerr << "\033[1;33m";
        cerr << input_word;
        cerr << "\033[0m";
        cerr << " disregarding this variable and continueing to next line";
        cerr << endl;
        getline(input_file, input_word);
    }
    check_parameters();
    
    input_file.close();
}


// Function calaculate the bond order paramter
// GSL
// Y_l,m(theta, phi) = gsl_sf_legendre_sphPlm(l, m, cos(theta))*(cos(m*phi)+i*sin(m*phi))
// BOOST
// Y_l,m(theta, phi) = boost::math::spherical_harmonic_r(bond_parameter_order_, i_m, theta, phi);
void BondOrderParameter::compute_BOP()
{
    if (!is_wrapped_) {
        wrap_coordinates();
    }
    
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
    // select the indexes of atom_types_
    size_t number_of_atoms = 0;
    vector < vector< unsigned int > > atom_types_indexes(atom_types_.size());
    
    determine_atom_indexes(atom_types_indexes, number_of_atoms);
    
    // check max_cutoff_length < box_length_/2.0
    double min_box_length = average_box_length_[0];
    for (size_t i_dimension = 1; i_dimension < dimension_; ++i_dimension) {
        min_box_length = (min_box_length < average_box_length_[i_dimension]) ? min_box_length : average_box_length_[i_dimension];
    }
    if (max_cutoff_length_ == 0.0) {
        max_cutoff_length_ = min_box_length/2.0;
    }
    if (max_cutoff_length_ > min_box_length/2.0) {
        cerr << "WARNING: max cutoff length is greater than half of the smallest box length" << endl;
        cerr << "       : Setting the max cutoff length to half of the smallest box length" << endl;
        max_cutoff_length_ = min_box_length/2.0;
    }
    
    // resize the bond_order_parameter_
    if (is_averaged_) {
        bond_order_parameter_.resize(number_of_time_points_, vector<double> (1, 0.0));
    }
    else {
        bond_order_parameter_.resize(number_of_time_points_, vector<double> (number_of_atoms, 0.0));
    }
    
    // declare vectors for the real and imaginary term from spherical harmonics
    vector< vector<double> > real_term;
    vector< vector<double> > imaginary_term;
    
    // resize the real and imaginary term
    if (add_bar_ && !is_averaged_) {
        real_term.resize(number_of_atoms, vector<double> (bond_parameter_order_ + 1, 0.0));
        imaginary_term.resize(number_of_atoms, vector<double> (bond_parameter_order_ + 1, 0.0));
    }
    else {
        real_term.resize(1, vector<double> (bond_parameter_order_ + 1, 0.0));
        imaginary_term.resize(1, vector<double> (bond_parameter_order_ + 1, 0.0));
    }

    size_t status = 0;
    cout << "\nComputing ..." << endl;
    
    // main loop
#pragma omp parallel for firstprivate(real_term, imaginary_term)
    for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {

        // reinitialize the vectors
        if (add_bar_ && !is_averaged_) {
            real_term         = vector< vector <double> > (number_of_system_atoms_, vector<double> (bond_parameter_order_ + 1, 0.0));
            imaginary_term    = vector< vector <double> > (number_of_system_atoms_, vector<double> (bond_parameter_order_ + 1, 0.0));
        }
        else {
            real_term[0]      = vector<double> (bond_parameter_order_ + 1, 0.0);
            imaginary_term[0] = vector<double> (bond_parameter_order_ + 1, 0.0);
        }
        
        size_t current_frame    = time_array_indexes_[time_point];
        size_t i_atom_number    = 0;
        size_t number_neighbors = 0;
        
        // loop over all the atoms
        for (size_t i_atom_type1 = 0; i_atom_type1 < atom_types_.size(); ++i_atom_type1) {
            for (size_t i_atom1 = 0; i_atom1 < atom_types_indexes[i_atom_type1].size(); ++i_atom1, ++i_atom_number) {
                size_t atom1_index = atom_types_indexes[i_atom_type1][i_atom1];

                if (!is_averaged_) {
                    number_neighbors = 0;
                    
                    if (!add_bar_) {
                        real_term[0]      = vector<double> (bond_parameter_order_ + 1, 0.0);
                        imaginary_term[0] = vector<double> (bond_parameter_order_ + 1, 0.0);
                    }
                }
                
                // loop over all the atoms again
                for (size_t i_atom_type2 = 0; i_atom_type2 < atom_types_.size(); ++i_atom_type2) {
                    for (size_t i_atom2 = 0; i_atom2 < atom_types_indexes[i_atom_type2].size(); ++i_atom2) {
                        size_t atom2_index = atom_types_indexes[i_atom_type2][i_atom2];
                        
                        if (atom1_index == atom2_index) {
                            continue;
                        }
                        
                        // compute the spherical harmonic terms
                        if (add_bar_ && !is_averaged_) {
                            compute_harmonic(atom1_index, atom2_index, current_frame, number_neighbors, real_term[atom1_index], imaginary_term[atom1_index]);
                        }
                        else {
                            compute_harmonic(atom1_index, atom2_index, current_frame, number_neighbors, real_term[0], imaginary_term[0]);
                        }
                    }
                }

                if (add_bar_ && !is_averaged_) {
                    if (number_neighbors != 0) {
                        for (size_t i_m = 0; i_m < bond_parameter_order_ + 1; ++i_m) {
                            real_term[atom1_index][i_m]      /= number_neighbors;
                            imaginary_term[atom1_index][i_m] /= number_neighbors;
                        }
                    }
                }
                else if (!is_averaged_) {
                    if (number_neighbors != 0) {
                        calculate_BOP(real_term[0], imaginary_term[0], number_neighbors, time_point, i_atom_number);
                    }
                }
            }
        }
        
        // need to reloop the atoms if we want q bar
        if (add_bar_ && !is_averaged_) {
            i_atom_number = 0;
            
            for (size_t i_atom_type1 = 0; i_atom_type1 < atom_types_.size(); ++i_atom_type1) {
                for (size_t i_atom1 = 0; i_atom1 < atom_types_indexes[i_atom_type1].size(); ++i_atom1, ++i_atom_number) {
                    size_t atom1_index = atom_types_indexes[i_atom_type1][i_atom1];

                    vector<double> real_term_bar      (bond_parameter_order_ + 1, 0.0);
                    vector<double> imaginary_term_bar (bond_parameter_order_ + 1, 0.0);
                    
                    number_neighbors = 0;
                    
                    for (size_t i_atom_type2 = 0; i_atom_type2 < atom_types_.size(); ++i_atom_type2) {
                        for (size_t i_atom2 = 0; i_atom2 < atom_types_indexes[i_atom_type2].size(); ++i_atom2) {
                            size_t atom2_index = atom_types_indexes[i_atom_type2][i_atom2];
                            
                            double distance = 0.0;
                            vector<double> delta_r(dimension_, 0.0);
                            
                            if (atom1_index != atom2_index) {
                                for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                    delta_r[i_dimension] = trajectory_[current_frame][atom2_index][i_dimension] - trajectory_[current_frame][atom1_index][i_dimension];
                                    delta_r[i_dimension] -= box_length_[current_frame][i_dimension] * round(delta_r[i_dimension]/box_length_[current_frame][i_dimension]);
                                    distance += delta_r[i_dimension]*delta_r[i_dimension];
                                }
                                distance = sqrt(distance);
                            }
                            
                            // average q over the nearest neighbors
                            if (distance < max_cutoff_length_) {
                                for (size_t i_m = 0; i_m < bond_parameter_order_ + 1; ++i_m) {
                                    real_term_bar[i_m]      += real_term[atom2_index][i_m];
                                    imaginary_term_bar[i_m] += imaginary_term[atom2_index][i_m];
                                }
                                
                                number_neighbors += 1;
                            }
                        }
                    }
                    
                    if (number_neighbors != 0) {
                        calculate_BOP(real_term_bar, imaginary_term_bar, number_neighbors, time_point, i_atom_number);
                    }
                }
            }
        }
        else if (is_averaged_) {
            if (number_neighbors != 0) {
                calculate_BOP(real_term[0], imaginary_term[0], number_neighbors, time_point, 0);
            }
        }
        
        if (is_run_mode_verbose_) {
#pragma omp critical
            {
                print_status(status);
            }
        }
    }

    cout << endl;
}


// method calculate and store bond_order_parameter
// q_l = q_l0*q_l0 + 2.0*sum_1:l q_li
void BondOrderParameter::calculate_BOP(vector<double> const & real_term, vector<double> const & imaginary_term,
                                       size_t const & number_neighbors, size_t const & time_point, size_t const & atom_index)
{
    bond_order_parameter_[time_point][atom_index] = (real_term[0]*real_term[0] + imaginary_term[0]*imaginary_term[0]);
    for (size_t i_m = 1; i_m < bond_parameter_order_ + 1; ++i_m) {
        bond_order_parameter_[time_point][atom_index] += 2.0*(real_term[i_m]*real_term[i_m] + imaginary_term[i_m]*imaginary_term[i_m]);
    }
    bond_order_parameter_[time_point][atom_index] = sqrt(bond_order_parameter_[time_point][atom_index]);
    bond_order_parameter_[time_point][atom_index] *= sqrt(4.0*M_PI/(2.0*bond_parameter_order_ + 1.0))/number_neighbors;
}


// method to perform the spherical harmonics
void BondOrderParameter::compute_harmonic(size_t const & atom1_index, size_t const & atom2_index, size_t const & current_frame,
                                          size_t & number_neighbors, vector<double> & real_term, vector<double> & imaginary_term)
{
    double distance = 0.0;
    vector<double> delta_r(dimension_,0);
    
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        delta_r[i_dimension] = trajectory_[current_frame][atom2_index][i_dimension] - trajectory_[current_frame][atom1_index][i_dimension];
        delta_r[i_dimension] -= box_length_[current_frame][i_dimension] * round(delta_r[i_dimension]/box_length_[current_frame][i_dimension]);
        distance += delta_r[i_dimension]*delta_r[i_dimension];
    }
    distance = sqrt(distance);
    
    if (distance < max_cutoff_length_) {
        number_neighbors += 1;
        
        double theta = acos(delta_r[2]/distance);
        double phi   = atan2(delta_r[1], delta_r[0]);
        
        // bond_order_degree is the term "m" in the spherical harmonics equation
        // i_m is the iterator over the vector, and bond_order_degree is the conversion to the order "m"
        for (size_t i_m = 0; i_m < bond_parameter_order_ + 1; ++i_m) {
        #ifdef GSL
            real_term[i_m]      += gsl_sf_legendre_sphPlm(bond_parameter_order_, i_m, cos(theta))*cos(i_m*phi);
            imaginary_term[i_m] += gsl_sf_legendre_sphPlm(bond_parameter_order_, i_m, cos(theta))*sin(i_m*phi);
        #elif BOOST
            real_term[i_m]      += boost::math::spherical_harmonic_r(bond_parameter_order_, i_m, theta, phi);
            imaginary_term[i_m] += boost::math::spherical_harmonic_i(bond_parameter_order_, i_m, theta, phi);
        #endif
        }
    }
}


void BondOrderParameter::write_BOP()
{
    // set the time step of the system if not defined
    if (trajectory_delta_time_ == 0.0) {
        cerr << endl;
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
    size_t number_of_atoms = 0;
    vector < vector< unsigned int > > atom_types_indexes(atom_types_.size());
    
    if (!is_averaged_) {
        // redict output so not printed to screen
        streambuf *standard = cout.rdbuf();
        stringstream temp_stream;
        cout.rdbuf(temp_stream.rdbuf());
        
        // select the indexes of atom_types_
        determine_atom_indexes(atom_types_indexes, number_of_atoms);
        
        // reset cout to print to screen
        cout.rdbuf(standard);
    }

    // open output file
    ofstream output_BOP_file(output_file_name_);
    
    // print the header of the file
    output_BOP_file    << "# Bonded Order Parameter for ";
    output_BOP_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_BOP_file << atom_types_[i_atom_type];
        output_BOP_file << " ";
    }
    output_BOP_file << "}";
    output_BOP_file << "in " << atom_group_ << ".\n";
    if (time_scale_type_ == "linear") {
        output_BOP_file << "# using " << time_scale_type_ << " scale\n";
    }
    else if (time_scale_type_ == "log") {
        output_BOP_file << "# using " << time_scale_type_ << " scale, logscale resulted in " << number_of_time_points_ - time_array_indexes_.size() << " repeated points ignored\n";
    }
    if (is_averaged_) {
        output_BOP_file << "# time | Q_" << bond_parameter_order_ << "\n";
        output_BOP_file << "#\n";
        output_BOP_file << setiosflags(ios::scientific) << setprecision(output_precision_);
        
        // print BOP
        for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
            output_BOP_file << time_array_indexes_[time_point]*trajectory_delta_time_;
            output_BOP_file << " ";
            output_BOP_file << bond_order_parameter_[time_point][0];
            output_BOP_file << "\n";
        }
    }
    else {
        if (add_bar_) {
            output_BOP_file << "# time | qbar_" << bond_parameter_order_ << "(atom)\n";
        }
        else {
            output_BOP_file << "# time | q_" << bond_parameter_order_ << "(atom)\n";
        }
        output_BOP_file << "# atom | ";
        for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
            for (size_t i_atom = 0; i_atom < atom_types_indexes[i_atom_type].size(); ++i_atom) {
                output_BOP_file << atom_types_indexes[i_atom_type][i_atom] << " ";
            }
        }
        output_BOP_file << endl;
        output_BOP_file << setiosflags(ios::scientific) << setprecision(output_precision_);
        
        //print BOP
        for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
            output_BOP_file << time_array_indexes_[time_point]*trajectory_delta_time_;
            for (size_t i_atom = 0; i_atom < number_of_atoms; ++i_atom) {
                output_BOP_file << " ";
                output_BOP_file << bond_order_parameter_[time_point][i_atom];
            }
            output_BOP_file << "\n";
        }
    }

    output_BOP_file.close();
}


void BondOrderParameter::check_parameters() throw()
{
    if (end_frame_ == 0) {
        cerr << "ERROR: We require more information to proceed, either end_frame or number_of_time_points\n";
        cerr << "       must be povided for us to continue.";
        cerr << endl;
        exit(1);
    }
    
    if (time_scale_type_ == "linear") {
        if (end_frame_ == 0) {
            end_frame_ = start_frame_ + number_of_time_points_*frame_interval_;
        }
        if (number_of_time_points_ == 0) {
            number_of_time_points_ = (end_frame_ - start_frame_)/frame_interval_;
        }
        if (number_of_time_points_*frame_interval_ > end_frame_ - start_frame_) {
            end_frame_ = start_frame_ + number_of_time_points_*frame_interval_;
            cerr << "WARNING: the number of frames required is greater then the number supplied\n";
            cerr << "       : setting end frame to minimum value allowed: ";
            cerr << end_frame_;
            cerr << endl;
        }
        if (frame_interval_ < 1) {
            cerr << "ERROR: frame_interval must be an integer greater than 0 for linear scale\n" << endl;
            exit(1);
        }
    }
    else if (time_scale_type_ == "log") {
        if (end_frame_ == 0) {
            end_frame_ = start_frame_ + pow(frame_interval_,number_of_time_points_);
        }
        if (number_of_time_points_ == 0) {
            number_of_time_points_ = static_cast<unsigned int>(log(end_frame_ - start_frame_)/log(frame_interval_));
            cout << number_of_time_points_ << endl;
        }
        if (static_cast<unsigned int>(pow(frame_interval_, number_of_time_points_) + 0.5) > end_frame_ - start_frame_) {
            end_frame_ = start_frame_ + static_cast<unsigned int>(pow(frame_interval_, number_of_time_points_) + 0.5);
            cerr << "WARNING: the number of frames required is greater then the number supplied\n";
            cerr << "         setting end frame to minimum value allowed: ";
            cerr << end_frame_;
            cerr << endl;
        }
        if (frame_interval_ <= 1) {
            cerr << "ERROR: frame_interval must be greater than 1.0 for logscale\n" << endl;
            exit(1);
        }
    }
    else {
        cerr << "ERROR: Illegal time scale specified. Must be one of (linear/log)\n" << endl;
        exit(1);
    }
    
    if (!is_wrapped_) {
        cerr << "WARNING: the trajectory provided is not wrapped\n";
        cerr << "       : We will wrapp it for you, but user discretion\n";
        cerr << "       : is advised";
        cerr << endl;
    }
    
    if (atom_types_.empty()) {
        cerr << "ERROR: No atom types provided, nothing will be computed\n";
        cerr << endl;
        exit(1);
    }
    
    if (bond_parameter_order_ < 1) {
        cerr << "ERROR: Bond order must be an integer greater than 0\n";
        cerr << endl;
        exit(1);
    }
    
    if (max_cutoff_length_ <= 0.0) {
        cerr << "ERROR: A postive non-zero cutoff length must be provided\n";
        cerr << "     : The recommended value is the distance to the first minimum of the pair distribution\n";
        cerr << endl;
        exit(1);
    }
    
    if (dimension_ != 3) {
        cerr << "ERROR: Dimension must be 3 for this quantity\n";
        cerr << endl;
        exit(1);
    }
    
    // check output file can be opened
    ofstream output_BOP_file(output_file_name_);
    if (!output_BOP_file) {
        cerr << "ERROR: Output file for mean squared displacement: "
        << "\033[1;25m"
        << output_file_name_
        << "\033[0m"
        << ", could not be opened"
        << endl;
        exit(1);
    }
    output_BOP_file.close();
}


void BondOrderParameter::compute_time_array()
{
    // first time point is always zero
    time_array_indexes_.push_back(0);
    
    double       total_time     = frame_interval_;
    unsigned int frame_previous = 0;
    
    // TODO switch to try catch
    for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {
        // check that we have not exceeded allowed frames
        assert(static_cast<unsigned int>(total_time) < end_frame_ - start_frame_ && "Error: Not eneough frames for calculation on log time scale");
        
        if (time_scale_type_ == "linear") {
            time_array_indexes_.push_back(static_cast<unsigned int>(total_time));
            total_time += frame_interval_;
        }
        else if (time_scale_type_ == "log") {
            if (static_cast<unsigned int>(total_time) != frame_previous) {
                time_array_indexes_.push_back(static_cast<unsigned int>(total_time));
                frame_previous = static_cast<unsigned int>(total_time);
            }
            total_time *= frame_interval_;
        }
    }
}


void BondOrderParameter::determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                                size_t & number_of_atoms)
{
    number_of_atoms = 0;
    
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        select_atoms(atom_types_indexes[i_atom_type], atom_types_[i_atom_type], atom_group_);
        number_of_atoms += atom_types_indexes[i_atom_type].size();
    }
}


void BondOrderParameter::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the bond order parameter is: ";
    cout << status * 100.0/time_array_indexes_.size();
    cout << " \%";
    cout << flush;
}