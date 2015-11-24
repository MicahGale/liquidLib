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
    method_of_k_generation_("moderate"),
    number_of_bins_(50),
    k_start_index_(0),
	number_of_frames_to_average_(0),
    atom_types_({}),
    scattering_lengths_({}),
    number_of_k_vectors_({})
{
}


StructureFactor::~StructureFactor()
{
}


void StructureFactor::read_command_inputs(int argc, char * argv[])
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
            is_run_mode_verbose_ = true;
            continue;
        }
        
        cerr << "\nERROR: Unrecognized flag '" << argv[input] << "' from command inputs.\n";
        exit(1);
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
			cerr << "     : compatible version of LiquidLib\n";
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
        if (input_word == "trajectory_data_type") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            trajectory_data_type_ = input_word;
            continue;
        }
        
        // Read scattering lengths provided by user
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
        if (input_word == "scattering_lengths") {
            getline(input_file, input_word);
            stringstream input_line(input_word);
            while (input_line >> input_word) {
                if (input_word[0] == '=') {
                    input_line >> input_word;
                }
                if (input_word[0] == '#') {
                    break;
                }
                scattering_lengths_.push_back(stod(input_word));
            }
            continue;
        }
        
        if (input_word == "method_of_k_generation") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            method_of_k_generation_ = input_word;
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
	
	cout << "\nInput file reading complete\n" << endl;
	input_file.close();
}


void StructureFactor::compute_S_k()
{
    if (!is_wrapped_) {
        wrap_coordinates();
    }
    
    // select the indexes of atom_types_
    size_t number_of_atoms;
    double average_scattering_length = 0.0;
    vector < vector< unsigned int > > atom_types_indexes(atom_types_.size());
    
    determine_atom_indexes(atom_types_indexes, average_scattering_length, number_of_atoms);
    
    // k resolution is determined by inverse box length
    double min_box_length = average_box_length_[0];
    for (size_t i_dimension = 1; i_dimension < dimension_; ++i_dimension) {
        min_box_length = (min_box_length < average_box_length_[i_dimension]) ? min_box_length : average_box_length_[i_dimension];
    }
    double delta_k = 2.0*M_PI/min_box_length;
    
    // allocate k_values_
    k_values_.resize(number_of_bins_, 0.0);
    number_of_k_vectors_.resize(number_of_bins_);
    vector< vector< vector < unsigned int > > > k_vectors(number_of_bins_);
    vector< double > normalization_factor(number_of_bins_, 0.0);
    
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        // Compute k_vectors
        k_values_[k_index] = delta_k * (k_index + k_start_index_);
        generate_k_vectors(k_index + k_start_index_, k_vectors[k_index]);
        number_of_k_vectors_[k_index] = k_vectors[k_index].size();
        
        // Compute normalization factor
        normalization_factor[k_index] = 1.0 / (number_of_frames_to_average_ * number_of_atoms * number_of_k_vectors_[k_index]);
        normalization_factor[k_index] /= ( average_scattering_length * average_scattering_length );
    }
    
    // allocate S_k_
    S_k_.resize(number_of_bins_, 0.0);
    
    size_t status = 0;
    cout << "Computing ..." << endl;
    
#pragma omp parallel for firstprivate(k_vectors)
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        for (size_t frame_number = 0; frame_number < number_of_frames_to_average_; ++frame_number) {
            for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_[k_index]; ++i_k_vector) {
                double sum_cos_term = 0.0;
                double sum_sin_term = 0.0;
                for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
                    for (size_t i_atom = 0; i_atom < atom_types_indexes[i_atom_type].size(); ++i_atom) {
                        size_t atom_index = atom_types_indexes[i_atom_type][i_atom];
                        double k_dot_r = 0.0;
                        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            k_dot_r += k_vectors[k_index][i_k_vector][i_dimension] * trajectory_[frame_number][atom_index][i_dimension];
                        }
                        sum_cos_term += scattering_lengths_[i_atom_type] * cos(k_dot_r);
                        sum_sin_term += scattering_lengths_[i_atom_type] * sin(k_dot_r);
                    }
                }
                S_k_[k_index] += (sum_cos_term * sum_cos_term + sum_sin_term * sum_sin_term);
            }
        }
        S_k_[k_index] *= normalization_factor[k_index];

        if (is_run_mode_verbose_) {
#pragma omp atomic
            print_status(status);
        }
    }
    cout << endl;
}


void StructureFactor::write_S_k()
{
	ofstream output_Sk_file(output_file_name_);
	
	output_Sk_file << "# Total structure factor for atom types: \n";
    output_Sk_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_Sk_file << atom_types_[i_atom_type];
        output_Sk_file << "(" << scattering_lengths_[i_atom_type] << ")";
        output_Sk_file << " ";
    }
    output_Sk_file << "}";
	output_Sk_file << "in " << atom_group_ << ".\n";
    output_Sk_file << "# k-values       S(k)  " << "\n";
    
    output_Sk_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
		output_Sk_file << k_values_[k_index];
        output_Sk_file << "       ";
        output_Sk_file << S_k_[k_index];
        output_Sk_file << "\n";
	}
    
	output_Sk_file.close();
}


void StructureFactor::generate_k_vectors(unsigned int const & mag_kvec_sqr, vector< vector< unsigned int > > & k_vectors)
{
    if (method_of_k_generation_ == "moderate") {
        unsigned int mag_k_vector = static_cast< unsigned int > (sqrt(mag_kvec_sqr));
        for (unsigned int i = 0; i <= mag_k_vector; ++i) {
            unsigned int j_max = static_cast< unsigned int > (sqrt(mag_kvec_sqr - i*i));
            for (unsigned int j = 0; j <= j_max; ++j) {
                unsigned int k_sqr = mag_kvec_sqr - i*i - j*j;
                unsigned int k = static_cast< unsigned int > (sqrt(k_sqr));
                // check if k is a perfect square integer
                if ( fabs(k_sqr - k*k) < 1e-16) {
                    vector< unsigned int > inds = {i, j, k};
                    k_vectors.push_back( inds );
                }
            }
        }
        return;
    }
    else if (method_of_k_generation_ == "fast") {
        unsigned int i_max = static_cast< unsigned int > (sqrt(mag_kvec_sqr/3.0));
        for (unsigned int i = 0; i <= i_max; ++i) {
            unsigned int j_max = static_cast< unsigned int > (sqrt((mag_kvec_sqr - i*i)/2.0));
            for (unsigned int j = i; j <= j_max; ++j) {
                unsigned int k_sqr = mag_kvec_sqr - i*i - j*j;
                unsigned int k = static_cast< unsigned int > (sqrt(k_sqr));
                // check if k is a perfect square
                if ( fabs(k_sqr - k*k) < 1e-16 ) {
                    // Add all permutations since i <= j <= k;
                    k_vectors.push_back( {i, j, k} );
                    if (i == j && j == k) {
                        continue;
                    }
                    else if ( i == j || j == k ) {
                        k_vectors.push_back( {j, k, i} );
                        k_vectors.push_back( {k, i, j} );
                    }
                    else {
                        // Write down all 6 permuations
                        k_vectors.push_back( {i, k, j} );
                        k_vectors.push_back( {j, i, k} );
                        k_vectors.push_back( {j, k, i} );
                        k_vectors.push_back( {k, i, j} );
                        k_vectors.push_back( {k, j, i} );
                    }
                }
            }
        }
    }
    else if (method_of_k_generation_ == "dimensional") {
        vector< unsigned int > k_vector(dimension_, 0);
        recurrsive_k_generation(mag_kvec_sqr, 0, k_vector, k_vectors);
    
        return;
    }
    // add other sampling methods that are also generic in different dimensions
}


void StructureFactor::recurrsive_k_generation(unsigned int const mag_kvec_sqr, unsigned int const dimension,
                                              vector< unsigned int > & k_vector, vector< vector< unsigned int > > & k_vectors)
{
    unsigned int sum_squares = 0;
    for (size_t i_dimension = 0; i_dimension < dimension; ++i_dimension) {
        sum_squares += k_vector[i_dimension]*k_vector[i_dimension];
    }
    unsigned int max_iterator = static_cast< unsigned int > (sqrt(mag_kvec_sqr - sum_squares));
    for (size_t k_mag = 0; k_mag <= max_iterator; ++k_mag) {
        k_vector[dimension] = k_mag;
        if (dimension != dimension_ - 1) {
            recurrsive_k_generation(mag_kvec_sqr, dimension + 1, k_vector, k_vectors);
        }
        else {
            unsigned int k_squared = 0;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                k_squared += k_vector[i_dimension]*k_vector[i_dimension];
            }
            if (k_squared == mag_kvec_sqr) {
                k_vectors.push_back(k_vector);
            }
        }
    }
}


// Function to check that all the parameters provided by the user for
// structure factor are usable.  Ensures that the code will
// not have with needing enough data points
void StructureFactor::check_parameters() throw()
{
    // Neither parameters exists
    if (number_of_frames_to_average_ == 0 && end_frame_ == 0) {
        cerr << "\n";
        cerr << "ERROR: Either 'number_of_frames_to_average' or 'end_frame' is not specified in input file.\n";
        cerr << "     : Not enough information to read trajectory" << endl;
        exit(1);
    }
    
    // Both parameters exist
    if (number_of_frames_to_average_ > 0 && end_frame_ > 0) {
        cerr << "\n";
        cerr << "WARNING: Both 'end_frame' and 'number_of_frames_to_average' values set in input file.\n";
        cerr << "       : 'number_of_frames_to_average' is used for reading necessary frames and for computation.";
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
    
    if (method_of_k_generation_ != "moderate" && method_of_k_generation_ != "fast" && method_of_k_generation_ != "dimensional") {
        cerr << "\n";
        cerr << "ERROR: Unrecognized sampling method for wavevector transfer k" << endl;
        cerr << "     : Optional methods are \"moderate\", \"fast\", or \"dimensional\"." << endl;
        exit(1);
    }
    
    if (dimension_ != 3 && method_of_k_generation_ != "dimensional") {
        cout << "WARNING: Only dimensional k_generation works for high dimensional trajectories\n";
        cout << "       : We will set k_generation_method to dimension\n";
        cout << endl;
        
        method_of_k_generation_ = "dimensional";
    }
    
    if (k_start_index_ < 1) {
        cerr << "WARNING: k_start_index cannot be non-positive\n";
        cerr << "       : We will reset the value to 1\n";
        cerr << endl;
        k_start_index_ = 1;
    }
    
    if (atom_types_.empty()) {
        cerr << "ERROR: No atom types provided, nothing will be computed\n";
        cerr << endl;
        exit(1);
    }
    
    // Check scattering lengths is correct
    if (scattering_lengths_.empty()) {
        cout << "Scattering lengths of atoms are not provided.\n";
        cout << "  This calculation will not be weighted by scattering length.\n";
        cout << endl;
    }
    else if (scattering_lengths_.size() != atom_types_.size() ) {
        cerr << "ERROR: Number of scattering lengths must match the number of atom types provided.\n";
        cerr << "     : Cannot proceed with computation\n";
        cerr << endl;
        exit(1);
    }
    
    if (scattering_lengths_.empty()) {
        scattering_lengths_ = vector< double > (atom_types_.size(), 1.0);
    }
    
    // check output file can be opened
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
    output_Sk_file.close();
}


void StructureFactor::determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                             double & average_scattering_length,
                                             size_t & number_of_atoms)
{
    number_of_atoms = 0;
    
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        select_atoms(atom_types_indexes[i_atom_type], atom_types_[i_atom_type], atom_group_);
        
        number_of_atoms += atom_types_indexes[i_atom_type].size();
        // Generate average_scattering_length
        if (scattering_lengths_.size() != 0) {
            average_scattering_length += atom_types_indexes[i_atom_type].size() * scattering_lengths_[i_atom_type];
        }
    }
    
    // normalize the average scattering length
    if (scattering_lengths_.size() == 0) {
        average_scattering_length = 1.0;
    }
    else {
        average_scattering_length /= number_of_atoms;
    }
}


void StructureFactor::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the pair distribution function is: ";
    cout << status * 100.0/number_of_bins_;
    cout << " \%";
    cout << flush;
}
