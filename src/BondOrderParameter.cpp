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

#include "Trajectory.hpp"

using namespace std;

BondOrderParameter::BondOrderParameter() :
    is_averaged_(false),
    input_file_name_("BOP.in"),
    output_file_name_("BOP.txt"),
    atom_group_("system"),
    time_scale_type_("linear"),
    number_of_time_points_(0),
    bond_parameter_order_(0),
    frame_interval_(1.0),
    max_cutoff_length_(0.0),
    atom_types_({})
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
// Y_l,m(theta, phi) = gsl_sf_legendre_sphPlm(l, m, cos(theta))*(cos(m*phi)+i*sin(m*phi))
void BondOrderParameter::compute_BOP()
{
#ifdef GSL
    if (!is_wrapped_) {
        wrap_coordinates();
    }
    
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
    // select the indexes of atom_types_
    size_t number_of_atoms;
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
    
    if (is_averaged_) {
        bond_order_parameter_.resize(number_of_time_points_, vector<double> (1, 0));
    }
    else {
        bond_order_parameter_.resize(number_of_time_points_, vector<double> (number_of_system_atoms_, 0));
    }

    size_t status = 0;
    cout << "\nComputing ..." << endl;
    

    cout << endl << gsl_sf_legendre_sphPlm(6, 6, .5)*cos(6*.5) <<" + i" << gsl_sf_legendre_sphPlm(6, 6, .5)*sin(6*.5) << endl;

    for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
        for (size_t i_atom_type1 = 0; i_atom_type1 < atom_types_.size(); ++i_atom_type1) {
            for (size_t i_atom1 = 0; i_atom1 < atom_types_indexes[i_atom_type1].size(); ++i_atom1) {
                for (size_t i_atom_type2 = 0; i_atom_type2 < atom_types_.size(); ++i_atom_type2) {
                    for (size_t i_atom2 = 0; i_atom2 < atom_types_indexes[i_atom_type2].size(); ++i_atom2) {
                        continue;
                    }
                }
            }
        }
        if (is_run_mode_verbose_) {
#pragma omp critical
            {
                print_status(status);
            }
        }
    }
#endif
}


void BondOrderParameter::write_BOP()
{
    if (trajectory_delta_time_ == 0.0) {
        cerr << endl;
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
    ofstream output_BOP_file(output_file_name_);
    
    output_BOP_file    << "#Bonded Order Parameter for ";
    output_BOP_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_BOP_file << atom_types_[i_atom_type];
        output_BOP_file << " ";
    }
    output_BOP_file << "}";
    output_BOP_file << "in " << atom_group_ << ".\n";
    output_BOP_file << "#using " << time_scale_type_ << " scale\n";
    if (is_averaged_) {
        output_BOP_file << "#time | Q_" << bond_parameter_order_ << "\n";
        output_BOP_file << setiosflags(ios::scientific) << setprecision(output_precision_);
        for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
            cout << time_point << endl;
            output_BOP_file << time_array_indexes_[time_point]*trajectory_delta_time_;
            output_BOP_file << " ";
            output_BOP_file << bond_order_parameter_[time_point][0];
            output_BOP_file << "\n";
        }
    }
    else {
        output_BOP_file << "#time | q_" << bond_parameter_order_ << "(atom)\n";
        output_BOP_file    << setiosflags(ios::scientific) << setprecision(output_precision_);
        for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
            output_BOP_file << time_array_indexes_[time_point]*trajectory_delta_time_;
            for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
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
        cerr << "ERROR: We require more information to proceed, either frameend or numberoftimepoints\n";
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
    time_array_indexes_.resize(number_of_time_points_);
    time_array_indexes_[0] = 0;
    
    double       total_time     = frame_interval_;
    unsigned int frame_previous = 0;
    // TODO switch to try catch
    for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {
        if (time_scale_type_ == "linear") {
            time_array_indexes_[time_point] = static_cast<unsigned int>(total_time);
            total_time += frame_interval_;
        }
        else {
            time_array_indexes_[time_point] = time_array_indexes_[time_point - 1];
            while (time_array_indexes_[time_point] == frame_previous) {
                time_array_indexes_[time_point] = static_cast<unsigned int>(total_time + 0.5);
                total_time *= frame_interval_;
            }
            frame_previous = time_array_indexes_[time_point];
        }
        
        assert(time_array_indexes_[time_point] < end_frame_ - start_frame_ && "Error: Not eneough frames for calculation on log time scale");
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
    cout << status * 100.0/number_of_time_points_;
    cout << " \%";
    cout << flush;
}