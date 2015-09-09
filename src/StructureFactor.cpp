//
//  SelfVanHoveFunction_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "StructureFactor.hpp"

#ifdef OMP
#include "omp.h"
#endif

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <iomanip>
#include <random>
#include <sstream>
#include <cstring>

#include "Trajectory.hpp"

using namespace std;

StructureFactor::StructureFactor() :
	input_file_name_("S_k.in"),
	output_file_name_("S_k.txt"),
	atom_group_("system"),
    method_of_k_sampling_("gaussian"),
    number_of_bins_(50),
    k_start_index_(0),
    number_of_k_vectors_(50),
	number_of_frames_to_average_(1)
{
}


StructureFactor::~StructureFactor()
{
}


void StructureFactor::read_command_inputs(int argc, char * argv[])
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


void StructureFactor::read_input_file()
{
	ifstream input_file(input_file_name_);
	
	if (!input_file) {
		cerr << "ERROR: Input file location, "
			 << "\033[1;25m"
			 << input_file_name_
			 << "\033[0m"
			 << ", does not exist, please check input"
			 << endl;
		exit(1);
	}
	
	cout << "Now reading the input file: "
		 << "\033[1;31m"
		 << input_file_name_
		 << "\033[0m"
		 << endl;
	
	string input_word;
	
	while (input_file >> input_word) {

		//check for comment
		if (input_word[0] == '#') {
			getline(input_file, input_word);
			continue;
		}
		
		//check for memeber bools
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
            if (output_file_name_ != "S_k.txt") {
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
			cerr << "ERROR: gro files cannot be used in non gromacs";
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
        
        // Read scattering lengths provided by user
		if (input_word == "scattering_lengths") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            input_line >> input_word;
            if (input_word[0] == '=') {
                input_line >> input_word;
            }
            while (input_word.find("#") == string::npos && !input_line.eof()) {
                user_atom_types_.push_back(input_word);
                input_line >> input_word;
                scattering_lengths_.push_back(stod(input_word));
                input_line >> input_word;
            }
            continue;
		}
        
        if (input_word == "method_of_k_sampling") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            method_of_k_sampling_ = input_word;
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
        if (input_word == "k_start_index") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            k_start_index_ = stoi(input_word);
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
		if (input_word == "number_of_k_vectors") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			number_of_k_vectors_ = stoi(input_word);
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
        if (input_word == "output_precision") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            output_precision_ = stod(input_word);
            continue;
        }
        
		//check for everything else

		cerr << input_word << endl;
		cerr << "ERROR: no matching input type for: ";
		cerr << "\033[1;33m";
		cerr << input_word;
		cerr << "\033[0m";
		cerr << " disregarding this variable and continueing to next line";
		cerr << endl;
		getline(input_file, input_word);
	}
	check_parameters();
	
	cout << "\nInput file reading complete" << endl;
	input_file.close();
}


void StructureFactor::compute_S_k()
{
    if (!is_wrapped_) {
        wrap_coordinates();
    }
    
    // k resolution is determined by inverse box length
    double min_box_length = average_box_length_[0];
    for (size_t i_dimension = 1; i_dimension < dimension_; ++i_dimension) {
        min_box_length = (min_box_length < average_box_length_[i_dimension]) ? min_box_length : average_box_length_[i_dimension];
    }
    double delta_k = 2.0*M_PI/min_box_length;
    
    // allocate k_values_ and S_k_
    k_values_.resize(number_of_bins_, 0.0);
    
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        k_values_[k_index] = delta_k * (k_index);
    }
    S_k_.resize(number_of_bins_, 0.0);
    
    // Compute normalization factor
    double average_scattering = 0.0;
    unsigned int total_number_of_atoms = 0;
    vector< vector< unsigned int > > atom_type_indexes(user_atom_types_.size());
    for (size_t i_atom_type = 0; i_atom_type < user_atom_types_.size(); ++i_atom_type) {
        select_atoms(atom_type_indexes[i_atom_type], user_atom_types_[i_atom_type], atom_group_);
        average_scattering += scattering_lengths_[i_atom_type] * atom_type_indexes[i_atom_type].size();
        total_number_of_atoms += atom_type_indexes[i_atom_type].size();
    }
    double normalization_factor = total_number_of_atoms / (number_of_frames_to_average_ * average_scattering * average_scattering * number_of_k_vectors_);
        
    int status = 0;
    
    // Perform frame averaging for S_k_
#pragma omp parallel
{
    vector< vector < double > > k_vectors(number_of_k_vectors_, vector< double >(dimension_, 0.0));

    for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
#pragma omp for
        for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
            generate_k_vectors(k_vectors, k_values_[k_index]);
            for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_; ++i_k_vector) {
                double sum_cos_term = 0.0;
                double sum_sin_term = 0.0;
                for (size_t i_atom_type = 0; i_atom_type < user_atom_types_.size(); ++i_atom_type) {
                    for (size_t i_atom = 0; i_atom < atom_type_indexes[i_atom_type].size(); ++i_atom) {
                        size_t atom_index = atom_type_indexes[i_atom_type][i_atom];
                        double k_dot_r = 0.0;
                        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            k_dot_r += k_vectors[i_k_vector][i_dimension] * trajectory_[frame_number][atom_index][i_dimension];
                        }
                        sum_cos_term += scattering_lengths_[i_atom_type] * cos(k_dot_r);
                        sum_sin_term += scattering_lengths_[i_atom_type] * sin(k_dot_r);
                    }
                }
                S_k_[k_index] += (sum_cos_term * sum_cos_term + sum_sin_term * sum_sin_term);
            }
        }
#pragma omp single
        {
            ++status;
            cout << "\rcurrent progress of calculating structure factor is: ";
            cout << static_cast< double >(status) * 100.0/number_of_frames_to_average_;
            cout << " \%";
            cout << flush;
        }
#pragma omp barrier
    }
}
    
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        S_k_[k_index] *= normalization_factor;
    }
}


void StructureFactor::write_S_k()
{
	ofstream output_Sk_file(output_file_name_);
	
	if (!output_Sk_file) {
		cerr << "ERROR: Output file for structure factor: "
			 << "\033[1;25m"
			 << output_file_name_
			 << "\033[0m"
			 << ", could not be opened"
			 << endl;
		exit(1);
	}
	
    output_Sk_file << setiosflags(ios::scientific) << setprecision(output_precision_);
	output_Sk_file << "# Total structure factor for atom types: ";
    for (size_t i_atom_type = 0; i_atom_type < user_atom_types_.size(); ++i_atom_type) {
        output_Sk_file << user_atom_types_[i_atom_type];
        output_Sk_file << "(" << scattering_lengths_[i_atom_type] << ")";
        output_Sk_file << " , ";
    }
	output_Sk_file << " in group: " << atom_group_ << ".\n";
    output_Sk_file << "# k-values       S(k)  " << "\n";
    
    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
		output_Sk_file << k_values_[k_index];
        output_Sk_file << "       ";
        output_Sk_file << S_k_[k_index];
        output_Sk_file << "\n";
	}
    
	cout << "Structure factor printing complete. \n" << endl;
	output_Sk_file.close();
}


void StructureFactor::generate_k_vectors(vector< vector< double > > & k_vectors, double const & k_absolute_value)
{
    if (method_of_k_sampling_ == "gaussian") {
        random_device seed;
        default_random_engine generator(seed());
        normal_distribution< double > random_number(0.0, 1.0); // (Mean, Stdev)
        vector< double > random_vector(dimension_, 0.0);
        for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_; ++i_k_vector) {
            double vector_length = 0.0;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                random_vector[i_dimension] = random_number(generator);
                vector_length += random_vector[i_dimension] * random_vector[i_dimension];
            }
            vector_length = sqrt(vector_length);
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                k_vectors[i_k_vector][i_dimension] = k_absolute_value * random_vector[i_dimension] / vector_length;
            }
        }
        return;
    }
    else if (method_of_k_sampling_ == "uniform") {
        random_device seed;
        default_random_engine generator(seed());
        uniform_real_distribution< double > random_number(0.0, 1.0); // (min, max)
        double phi, theta, x, y, z;
        for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_; ++i_k_vector) {
            phi = 2.0*M_PI*random_number(generator);
            theta = acos(1.0 - 2.0*random_number(generator));
            x = sin(theta)*cos(phi);
            y = sin(theta)*sin(phi);
            z = cos(theta);
            k_vectors[i_k_vector][0] = x * k_absolute_value;
            k_vectors[i_k_vector][1] = y * k_absolute_value;
            k_vectors[i_k_vector][2] = z * k_absolute_value;
        }
        return;
    }
    // add other sampling methods that are also generic in different dimensions
}


// Function to check that all the parameters provided by the user for
// structure factor are usable.  Ensures that the code will
// not have with needing enough data points
void StructureFactor::check_parameters() throw()
{
    if (end_frame_ == 0) {
        end_frame_ = number_of_frames_to_average_;
    }
    
    if (number_of_frames_to_average_ == 1 && end_frame_ == 0) {
        cerr << "WARNING: frameend and number_of_frames_to_average\n";
        cerr << "         values not set in input file.\n";
        cerr << "         1 frame average will be done.";
        cerr << endl;
        end_frame_ = 1;
    }
    
    if (number_of_frames_to_average_ > 1 && end_frame_ > 0) {
        cerr << "WARNING: Both 'frameend' and 'number_of_frames_to_average'\n";
        cerr << "         values set in input file.\n";
        cerr << "         'number_of_frames_to_average' will be used.";
        cerr << endl;
        end_frame_ = number_of_frames_to_average_ - 1;
    }
    
    if (number_of_frames_to_average_ == 1 && end_frame_ > 0) {
        cerr << "WARNING: number_of_frames_to_average not defined\n";
        cerr << "         in input file.\n";
        cerr << "         Value will not be set to 'frameend' - 'framestart'.";
        cerr << endl;
        number_of_frames_to_average_ = end_frame_ - start_frame_;
    }
    
    if (user_atom_types_.empty() || scattering_lengths_.empty()) {
        cerr << "ERROR: No atom types and scattering lengths specified.\n";
        cerr << "       Computation cannot proceed.";
        cerr << endl;
        exit(1);
    }
    
    if (method_of_k_sampling_ != "gaussian" && method_of_k_sampling_ != "uniform") {
        cerr << "ERROR: Unrecognized sampling method for wavevector transfer k" << endl;
        cerr << "     : Optional methods are \"analytical\", \"gaussian\", or \"uniform\"." << endl;
        exit(1);
    }
    
    if (method_of_k_sampling_ == "uniform") {
        if (dimension_ > 3) {
            cerr << "ERROR: uniform sampling can only be used for 2 or 3 dimensions";
            cerr << endl;
            exit(1);
        }
    }
}
