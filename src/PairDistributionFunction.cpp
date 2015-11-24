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
#include <sstream>

#include "Trajectory.hpp"

using namespace std;

PairDistributionFunction::PairDistributionFunction() :
	input_file_name_("g_r.in"),
	output_file_name_("g_r.txt"),
    atom_types1_({}),
    scattering_lengths1_({}),
    atom_group1_("system"),
    atom_types2_({}),
    scattering_lengths2_({}),
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

        
        //check for memeber bools
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
		if (input_word == "atom_types1") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                atom_types1_.push_back(input_word);
            }
            continue;
		}
        if (input_word == "scattering_lengths1") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                scattering_lengths1_.push_back(stod(input_word));
            }
            continue;
        }
        if (input_word == "atom_group1") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            atom_group1_ = input_word;
            continue;
        }
		if (input_word == "atom_types2") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                atom_types2_.push_back(input_word);
            }
            continue;
		}
        if (input_word == "scattering_lengths2") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                scattering_lengths2_.push_back(stod(input_word));
            }
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
    // check if atom_types are the same for decreased computing time
    bool same_atom_types = check_atom_types_equivalent();
    
	// select the indexes of atom_types1_
    size_t number_of_atoms1 = 0;
    double average_scattering_length1 = 0.0;
	vector < vector< unsigned int > > atom_types1_indexes(atom_types1_.size());

    determine_atom_indexes(atom_types1_, scattering_lengths1_, atom_group1_, atom_types1_indexes, average_scattering_length1, number_of_atoms1);
    
    // select the indexes of atom_types2_
    size_t number_of_atoms2 = 0;
    double average_scattering_length2 = 0.0;
    vector < vector< unsigned int > > atom_types2_indexes(atom_types2_.size());
    
    if (same_atom_types) {
        atom_types2_indexes        = atom_types1_indexes;
        average_scattering_length2 = average_scattering_length1;
        number_of_atoms2           = number_of_atoms1;
    }
    else {
        determine_atom_indexes(atom_types2_, scattering_lengths2_, atom_group2_, atom_types2_indexes, average_scattering_length2, number_of_atoms2);
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
	
    size_t status = 0;
    cout << "Computing ..." << endl;
    
	if (same_atom_types) {
#pragma omp parallel for
		for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
            for (size_t i_atom_type1 = 0; i_atom_type1 < atom_types1_.size(); ++i_atom_type1) {
                for (size_t i_atom1 = 0; i_atom1 < atom_types1_indexes[i_atom_type1].size(); ++i_atom1) {
                    for (size_t i_atom_type2 = 0; i_atom_type2 < atom_types2_.size(); ++i_atom_type2) {
                        for (size_t i_atom2 = i_atom1 + 1; i_atom2 < atom_types2_indexes[i_atom_type2].size(); ++i_atom2) {
                            size_t atom1_index = atom_types1_indexes[i_atom_type1][i_atom1];
                            size_t atom2_index = atom_types2_indexes[i_atom_type2][i_atom2];
                            
                            // only some trajectories provide molecule id's, this checks if molecule_id was created
                            if (!molecule_id_.empty()) {
                                if (molecule_id_[atom1_index] == molecule_id_[atom2_index]) {
                                    continue;
                                }
                            }
                            histogram_g_r(frame_number, atom1_index, atom2_index, i_atom_type1, i_atom_type2, delta_r);
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
	}
	else {
#pragma omp parallel for
        for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
            for (size_t i_atom_type1 = 0; i_atom_type1 < atom_types1_.size(); ++i_atom_type1) {
                for (size_t i_atom1 = 0; i_atom1 < atom_types1_indexes[i_atom_type1].size(); ++i_atom1) {
                    for (size_t i_atom_type2 = 0; i_atom_type2 < atom_types2_.size(); ++i_atom_type2) {
                        for (size_t i_atom2 = 0; i_atom2 < atom_types2_indexes[i_atom_type2].size(); ++i_atom2) {
                            size_t atom1_index = atom_types1_indexes[i_atom_type1][i_atom1];
                            size_t atom2_index = atom_types2_indexes[i_atom_type2][i_atom2];
                            
                            // only some trajectories provide molecule id's, this checks if molecule_id was created
                            if (!molecule_id_.empty()) {
                                if (molecule_id_[atom1_index] == molecule_id_[atom2_index]) {
                                    continue;
                                }
                            }
                            histogram_g_r(frame_number, atom1_index, atom2_index, i_atom_type1, i_atom_type2, delta_r);
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
	}
    cout << endl;
	
	// do normalization
	double average_volume = 1.0;
	for (vector< double >::iterator i_boxlength = average_box_length_.begin(); i_boxlength != average_box_length_.end(); ++i_boxlength) {
		average_volume *= *i_boxlength;
	}
	
    double density_of_atom_type2 = number_of_atoms2 / average_volume;
    int scaling_factor = 1;
    if (same_atom_types) {
        scaling_factor = 2;
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
		double normalization_factor = 1.0 / (density_of_atom_type2 * volume_of_shell * number_of_atoms1 * number_of_frames_to_average_);
        normalization_factor /= (average_scattering_length1 * average_scattering_length2);
		g_r_[i_bin] *= normalization_factor * scaling_factor;
    }
}


void PairDistributionFunction::write_g_r()
{
	ofstream output_gr_file(output_file_name_);
	
    output_gr_file << setiosflags(ios::scientific) << setprecision(output_precision_);
	output_gr_file << "# Pair Distribution Function between atom types:" << endl;
    output_gr_file << "# 1st: ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types1_.size(); ++i_atom_type) {
        output_gr_file << atom_types1_[i_atom_type];
        output_gr_file << "(" << scattering_lengths1_[i_atom_type] << ")";
        output_gr_file << " ";
    }
	output_gr_file << "in " << atom_group1_ << endl;
	output_gr_file << "# 2nd: ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types2_.size(); ++i_atom_type) {
        output_gr_file << atom_types2_[i_atom_type];
        output_gr_file << "(" << scattering_lengths2_[i_atom_type] << ")";
        output_gr_file << " ";
    }
    output_gr_file << "in " << atom_group2_ << endl;
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
    // Neither parameters exists
    if (number_of_frames_to_average_ == 0 && end_frame_ == 0) {
        cerr << "\n";
        cerr << "ERROR: Either 'number_of_frames_to_average' or 'end_frame' is not specified in input file.\n";
        cerr << "       Not enough information to read trajectory\n" << endl;
        exit(1);
    }
    
    // Both parameters exist
    if (number_of_frames_to_average_ > 0 && end_frame_ > 0) {
        cerr << "\n";
        cerr << "WARNING: Both 'end_frame' and 'number_of_frames_to_average' values set in input file.\n";
        cerr << "         'number_of_frames_to_average' is used for reading necessary frames and for computation.\n";
        cerr << endl;
        end_frame_ = start_frame_ + number_of_frames_to_average_;
    }
    
    // number_of_frames_to_average doesn't exist
    if (number_of_frames_to_average_ == 0) {
        number_of_frames_to_average_ = end_frame_ - start_frame_;
    }
    
    // end_frame doesn't exist
    if (end_frame_ == 0) {
        end_frame_ = start_frame_ + number_of_frames_to_average_;
    }
    
    // atom_types1 must contain atleast one atom type
    if (atom_types1_.empty()) {
        cerr << "\nERROR: 'atom_type1' not specifed. Must supply the type of atoms to compute g_r from.\n" << endl;
        exit(1);
    }
    
    // atom_types2 must contain atleast one atom type
    if (atom_types2_.empty()) {
        cerr << "\nERROR: 'atom_type2' not specifed. Must supply the second type of atoms to compute with.\n" << endl;
        exit(1);
    }
    
    if (scattering_lengths1_.empty()) {
        cout << "Scattering lengths of atoms are not provided.\n";
        cout << "  This calculation will not be weighted by scattering length.\n";
        cout << endl;
    }
    else if (scattering_lengths1_.size() != atom_types1_.size() ) {
        cerr << "ERROR: Number of scattering lengths must match the number of atom types for atom_types1 provided.\n";
        cerr << "     : Cannot proceed with computation\n";
        cerr << endl;
        exit(1);
    }
    
    if (!scattering_lengths2_.empty() && (scattering_lengths2_.size() != atom_types2_.size()) ) {
        cerr << "ERROR: Number of scattering lengths must match the number of atom types for atom_types2 provided.\n";
        cerr << "     : Cannot proceed with computation\n";
        cerr << endl;
        exit(1);
    }
    
    if (!scattering_lengths1_.empty() == scattering_lengths2_.empty()) {
        cerr << "ERROR: if scattering lengths are provided for atom_types1 or atom_types2,\n";
        cerr << "     : scattering lengths must be provided for the other as well\n";
        cerr << endl;
        exit(1);
    }
    
    if (atom_group2_ == "") {
        cerr << "\nWARNING: 'atom_group2' not specifed. It is set the same as 'atom_group'.\n" << endl;
        atom_group2_ = atom_group1_;
    }
    
    if (scattering_lengths1_.empty()) {
        scattering_lengths1_ = vector< double > (atom_types1_.size(), 1.0);
        scattering_lengths2_ = vector< double > (atom_types2_.size(), 1.0);
    }
    
    // check that output file can be opened
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
    output_gr_file.close();
}


void PairDistributionFunction::histogram_g_r(size_t const & frame_number,
                                             size_t const & atom1_index, size_t const & atom2_index,
                                             size_t const & atom1_type, size_t const & atom2_type,
                                             double const & delta_r)
{
    // determine the real distrance between two atoms
    double distrance_of_two_atoms = 0.0;
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        double delta_x = trajectory_[frame_number][atom1_index][i_dimension] - trajectory_[frame_number][atom2_index][i_dimension];
        delta_x -= box_length_[frame_number][i_dimension] * round(delta_x/box_length_[frame_number][i_dimension]);
        distrance_of_two_atoms += delta_x * delta_x;
    }
    distrance_of_two_atoms = sqrt(distrance_of_two_atoms);
    
    unsigned int bin = round(distrance_of_two_atoms/delta_r);
    if (bin < number_of_bins_) {
#pragma omp atomic
        g_r_[bin] += ( scattering_lengths1_[atom1_type] * scattering_lengths2_[atom2_type] );
    }
}


void PairDistributionFunction::determine_atom_indexes(vector < string > const & atom_types,
                                                      vector < double > const & scattering_lengths,
                                                      string            const & atom_group,
                                                      vector < vector < unsigned int > > & atom_types_indexes,
                                                      double & average_scattering_length,
                                                      size_t & number_of_atoms)
{
    number_of_atoms = 0;
    
    for (size_t i_atom_type = 0; i_atom_type < atom_types.size(); ++i_atom_type) {
        select_atoms(atom_types_indexes[i_atom_type], atom_types[i_atom_type], atom_group);
        
        number_of_atoms += atom_types_indexes[i_atom_type].size();
        // Generate average_scattering_length
        if (scattering_lengths.size() != 0) {
            average_scattering_length += atom_types_indexes[i_atom_type].size() * scattering_lengths[i_atom_type];
        }
    }
    
    // normalize the average scattering length
    if (scattering_lengths.size() == 0) {
        average_scattering_length = 1.0;
    }
    else {
        average_scattering_length /= number_of_atoms;
    }
}


bool PairDistributionFunction::check_atom_types_equivalent()
{
    if (atom_group1_ != atom_group2_) {
        return false;
    }
    
    if (atom_types1_.size() != atom_types2_.size()) {
        return false;
    }
    
    for (size_t i_atom_type = 0; i_atom_type < atom_types1_.size(); ++i_atom_type) {
        if (find(atom_types2_.begin(), atom_types2_.end(), atom_types1_[i_atom_type]) == atom_types2_.end()) {
            return false;
        }
    }
    
    return true;
}


void PairDistributionFunction::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the pair distribution function is: ";
    cout << status * 100.0/number_of_frames_to_average_;
    cout << " \%";
    cout << flush;
}