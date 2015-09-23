//
// PairDistributionFunction.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "PairDistributionFunction.hpp"

#ifdef OMP
#include "omp.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cstring>
#include <algorithm>

#include "Trajectory.hpp"

using namespace std;

PairDistributionFunction::PairDistributionFunction() :
	input_file_name_("g_r.in"),
	output_file_name_("g_r.txt"),
	atom_type_("all"),
	atom_type2_(""),
	atom_group_("system"),
	atom_group2_(""),
	number_of_bins_(200),
	number_of_frames_to_average_(0),
	max_cutoff_length_(0.0)
{
}


PairDistributionFunction::~PairDistributionFunction()
{
}


void PairDistributionFunction::read_command_inputs(int argc, char * argv[])
{
    for (int input = 0; input < argc; ++input) {
        if (strcmp(argv[input], "-i") == 0) {
            input_file_name_ = argv[++input];
        }
        if (strcmp(argv[input], "-o") == 0) {
            output_file_name_ = argv[++input];
        }
        if (strcmp(argv[input], "-t") == 0) {
            trajectory_file_name_ = argv[++input];
        }
    }
}


void PairDistributionFunction::read_input_file()
{
	ifstream input_file(input_file_name_);
	
	if (!input_file) {
		cerr << "ERROR: input file location, ";
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
		
		//check if equal to member strings
        if (input_word == "output_file_name") {
            if (output_file_name_ != "g_r.txt") {
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
		if (input_word == "atom_type") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			atom_type_ = input_word;
			continue;
		}
		if (input_word == "atom_type2") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			atom_type2_ = input_word;
			continue;
		}
		if (input_word == "atom_group") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			atom_group_ = input_word;
			continue;
		}
		if (input_word == "atom_group2") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			atom_group2_ = input_word;
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
		if (input_word == "number_of_bins") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			number_of_bins_ = stoi(input_word);
			continue;
		}
		if (input_word == "number_of_frames_to_average") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			number_of_frames_to_average_ = stoi(input_word);
			continue;
		}

		//check if equal to member doubles
		if (input_word == "max_cutoff_length") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			max_cutoff_length_ = stod(input_word);
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


void PairDistributionFunction::compute_g_r()
{
	// select the indexes of atom_type1 and atom_type2
	vector< unsigned int > atom_type1_indexes;
	vector< unsigned int > atom_type2_indexes;
	select_atoms(atom_type1_indexes, atom_type_, atom_group_);
	if (atom_type_ != atom_type2_ || atom_group_ != atom_group2_) {
		select_atoms(atom_type2_indexes, atom_type2_, atom_group2_);
	}

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
    
	double delta_r = max_cutoff_length_/number_of_bins_;
	
	g_r_.resize(number_of_bins_, 0.0);
	
    int status = 0;
    
	if (atom_type_ == atom_type2_ && atom_group_ == atom_group2_) {
#pragma omp parallel for
		for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
			for (vector< unsigned int >::iterator i_atom1 = atom_type1_indexes.begin(); i_atom1 != atom_type1_indexes.end(); ++i_atom1) {
				for (vector< unsigned int >::iterator i_atom2 = i_atom1 + 1; i_atom2 != atom_type1_indexes.end(); ++i_atom2) {
                    // only some trajectories provide molecule id's, this checks if molecule_id was created
                    if (!molecule_id_.empty()) {
                        if (molecule_id_[*i_atom1] == molecule_id_[*i_atom2]) {
                            continue;
                        }
                    }
                    histogram_g_r(frame_number, *i_atom1, *i_atom2, delta_r);
				}
			}
#pragma omp critical
{
            ++status;
            cout << "\rcurrent progress of calculating the pair distribution function is: ";
            cout << status * 100.0/number_of_frames_to_average_;
            cout << " \%";
            cout << flush;
}
		}
	}
	else {
#pragma omp parallel for
		for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
			for (vector< unsigned int >::iterator i_atom1 = atom_type1_indexes.begin(); i_atom1 != atom_type1_indexes.end(); ++i_atom1) {
				for (vector< unsigned int >::iterator i_atom2 = atom_type2_indexes.begin(); i_atom2 != atom_type2_indexes.end(); ++i_atom2) {
                    // only some trajectories provide molecule id's, this checks if molecule_id was created
                    if (!molecule_id_.empty()) {
                        if (molecule_id_[*i_atom1] == molecule_id_[*i_atom2]) {
                            continue;
                        }
                    }
                    histogram_g_r(frame_number, *i_atom1, *i_atom2, delta_r);
				}
			}
#pragma omp critical
{
            ++status;
            cout << "\rcurrent progress of calculating the pair distribution function is: ";
            cout << status * 100.0/number_of_frames_to_average_;
            cout << " \%";
            cout << flush;
}
		}
	}
    cout << endl;
	
	// do normalization
	double average_volume = 1.0;
	for (vector< double >::iterator i_boxlength = average_box_length_.begin(); i_boxlength != average_box_length_.end(); ++i_boxlength) {
		average_volume *= *i_boxlength;
	}
	
    double density_of_atom_type2 = atom_type2_indexes.size() / average_volume;
    int scaling_factor = 1;
    if (atom_type_ == atom_type2_ && atom_group_ == atom_group2_) {
        scaling_factor = 2;
        density_of_atom_type2 = atom_type1_indexes.size() / average_volume;
    }
    
	r_values_.resize(number_of_bins_, 0.0);
	double dimension_scaling_factor = pow(M_PI, dimension_/2.0) / tgamma(1 + dimension_/2.0);
    
	for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
		r_values_[i_bin] = delta_r * i_bin;
		
		// compute shell volume in a general dimension
		double volume_of_outer_sphere = pow(r_values_[i_bin] + delta_r/2.0, dimension_) * dimension_scaling_factor;
		double volume_of_inner_sphere = 0.0;
		if (i_bin != 0) {
		    volume_of_inner_sphere = pow(r_values_[i_bin] - delta_r/2.0, dimension_) * dimension_scaling_factor;
		}
		double volume_of_shell = volume_of_outer_sphere - volume_of_inner_sphere;
		double normalization_factor = 1.0 / (density_of_atom_type2 * volume_of_shell * atom_type1_indexes.size() * number_of_frames_to_average_);
		g_r_[i_bin] *= normalization_factor * scaling_factor;
	}
}


void PairDistributionFunction::write_g_r()
{
	ofstream output_gr_file(output_file_name_);
	
	if (!output_gr_file) {
		cerr << "Output file for Pair distribution function: ";
		cerr << "\033[1;25m";
        cerr << output_file_name_;
		cerr << "\033[0m";
		cerr << ", could not be opened";
		cerr << endl;
		exit(1);
	}
	
    output_gr_file << setiosflags(ios::scientific) << setprecision(output_precision_);
	output_gr_file << "# Pair Distribution Function between two atom types:" << endl;
	output_gr_file << "# 1st: " << atom_type_ << " in " << atom_group_ << endl;
	output_gr_file << "# 2nd: " << atom_type2_ << " in " << atom_group2_ << endl;
	output_gr_file << "# r                 g(r)" << endl;
	
	for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
        output_gr_file << r_values_[i_bin];
        output_gr_file << "        ";
        output_gr_file << g_r_[i_bin];
        output_gr_file << endl;
	}
	
	output_gr_file.close();
}


void PairDistributionFunction::check_parameters() throw()
{
    // check the number of frames to average
    if (number_of_frames_to_average_ > end_frame_ - start_frame_) {
        cerr << "WARNING: Not enough frames to average" << endl;
        cerr << "       : Reset the number of frames to average to the total number of frames available";
        number_of_frames_to_average_ = end_frame_ - start_frame_;
    }
    
    if (number_of_frames_to_average_ == 0) {
        number_of_frames_to_average_ = end_frame_ - start_frame_;
    }
    
    if (atom_type2_ == "") {
        atom_type2_ = atom_type_;
    }
    
    if (atom_group2_ == "") {
        atom_group2_ = atom_group_;
    }
}


void PairDistributionFunction::histogram_g_r(size_t const & frame_number, size_t const & i_atom1, size_t const & i_atom2, double const & delta_r)
{
    // determine the real distrance between two atoms
    double distrance_of_two_atoms = 0.0;
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        double delta_x = trajectory_[frame_number][i_atom1][i_dimension] - trajectory_[frame_number][i_atom2][i_dimension];
        delta_x -= box_length_[frame_number][i_dimension] * round(delta_x/box_length_[frame_number][i_dimension]);
        distrance_of_two_atoms += delta_x * delta_x;
    }
    distrance_of_two_atoms = sqrt(distrance_of_two_atoms);
    
    unsigned int bin = round(distrance_of_two_atoms/delta_r);
    if (bin < number_of_bins_) {
#pragma omp atomic
        g_r_[bin] += 1.0;
    }
}