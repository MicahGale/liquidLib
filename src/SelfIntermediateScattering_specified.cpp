//
//  SelfIntermediateScattering.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "SelfIntermediateScattering.hpp"

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
#include <cstring>
#include <sstream>

#include "Trajectory.hpp"

using namespace std;

SelfIntermediateScattering::SelfIntermediateScattering() :
	input_file_name_("Fs_kt.in"),
	output_file_name_("Fs_kt.txt"),
	atom_group_("system"),
	time_scale_type_("linear"),
    number_of_bins_(50),
    k_start_index_(0),
	number_of_time_points_(0),
	number_of_frames_to_average_(1),
	frame_interval_(1.0),
    delta_k_(0),
    atom_types_({}),
    scattering_lengths_({}),
    number_of_k_vectors_({})
{
}


SelfIntermediateScattering::~SelfIntermediateScattering()
{
}


void SelfIntermediateScattering::read_command_inputs(int argc, char * argv[])
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


void SelfIntermediateScattering::read_input_file()
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
            if (output_file_name_ != "Fs_kt.txt") {
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
        if (input_word == "delta_k") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            delta_k_ = stod(input_word);
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
		if (input_word == "number_of_frames_to_average") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			number_of_frames_to_average_ = stoi(input_word);
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


void SelfIntermediateScattering::compute_Fs_kt()		// Fs_kt is the Self intermediate scattering function
{
	if (is_wrapped_) {
		unwrap_coordinates();
	}
	
    // Form a array of time index values for a given type of timescale computation
	////////////////////////////////////////////////////////////////////////
	// is fine
	////////////////////////////////////////////////////////////////////////
    compute_time_array();
    
    // select the indexes of atom_types_
	////////////////////////////////////////////////////////////////////////
	// is fine
	////////////////////////////////////////////////////////////////////////
    size_t number_of_atoms;
    double average_squared_scattering_length = 0.0;
    vector < vector< unsigned int > > atom_types_indexes(atom_types_.size());
    
    determine_atom_indexes(atom_types_indexes, average_squared_scattering_length, number_of_atoms);

    // k resolution is determined by inverse box length
	////////////////////////////////////////////////////////////////////////
	// kill kill kill
	////////////////////////////////////////////////////////////////////////
    double min_box_length = average_box_length_[0];
    for (size_t i_dimension = 1; i_dimension < dimension_; ++i_dimension) {
        min_box_length = (min_box_length < average_box_length_[i_dimension]) ? min_box_length : average_box_length_[i_dimension];
    }
    double k_min = 2.0*M_PI/min_box_length;
    delta_k_ = (delta_k_ < k_min ) ? k_min : delta_k_;

    // allocate k_values_
	////////////////////////////////////////////////////////////////////////
	// kill kill kill
	// 
	// But also, you'll need to allocate space for the array eventually
	////////////////////////////////////////////////////////////////////////
    k_values_.resize(number_of_bins_, 0.0);
    number_of_k_vectors_.resize(number_of_bins_);
    vector< vector< vector < unsigned int > > > k_vectors(number_of_bins_);
    vector< double > normalization_factor(number_of_bins_, 0.0);
    
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
		////////////////////////////////////////////////////////////////////////
		// This section both generates the list of k-values and computes their k-vectors
		// So... we only need half of that
		// You probably want to just want to completely remake this for-loop
		// Pseudo-code:
		//
		// Allocate array space for k vectors, name it k_vectors to match remainder of code
		// FOR j in input_list_of_k_values:
		//      k_squared = integer((input_list_of_k_values[j])^2)
		//      compute_k_vectors(k_squared,k_vectors[k_index])
		//
		//      // Everything else here is (mostly) the same as in original loop
		//      // Skip the "k_values_[k_index] = blah blah blah" line
		//      number_of_k_vectors_[k_index] = k_vectors[k_index].size()
		//      normalization_factor[k_index] = 1.0 / (number_of_frames_to_average_ * number_of_atoms * number_of_k_vectors_[k_index]);
        //      normalization_factor[k_index] /= ( average_squared_scattering_length );
		//      // Copy pasta the rest of the cout lines to print the k-values to file
		//      
		////////////////////////////////////////////////////////////////////////
        
        // Compare delta_k_ with k_min allowed to select the larger of the two for first k-point
        double k_value_temp = delta_k_ * (k_index + k_start_index_);
        
        // Compute the integer value based on user input k_value
        unsigned int k_squared = static_cast<unsigned int> (k_value_temp * k_value_temp * min_box_length * min_box_length / (4.0 * M_PI * M_PI));
        
        // Compute the avaible k_vectors
        compute_k_vectors( k_squared, k_vectors[k_index] );
        number_of_k_vectors_[k_index] = k_vectors[k_index].size();
        
        // Corrected k_value on the grid with the new mag_kvec_sqr
        k_values_[k_index] = 2.0 * M_PI * sqrt(1.0 * k_squared) / min_box_length;
        
        // Compute normalization factor
        normalization_factor[k_index] = 1.0 / (number_of_frames_to_average_ * number_of_atoms * number_of_k_vectors_[k_index]);
        normalization_factor[k_index] /= ( average_squared_scattering_length );
        
        // Print num of vectors and values
        cout << "Corrected kvalue for : " << k_value_temp << " is : " << k_values_[k_index] << endl;
        cout << "Number of k-vectors : " << number_of_k_vectors_[k_index] << endl;
        if (is_run_mode_verbose_) {
            for (size_t i_vec = 0; i_vec < number_of_k_vectors_[k_index]; ++i_vec) {
                for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                    cout << k_vectors[k_index][i_vec][i_dimension] << ",";
                }
                cout << endl;
            }
        }
    }
    
    // initialize Fs_kt_
    Fs_kt_.resize(number_of_bins_, vector< double >(time_array_indexes_.size(), 0.0));   // initialize with value 1.0 since Fs(k, t = 0) = 1.0
    
    size_t status = 0;
    cout << "Computing ..." << endl;
    
#pragma omp parallel for firstprivate(k_vectors)
    for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
        for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {

            for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
                size_t current_frame = initial_frame + time_array_indexes_[time_point];

                for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_[k_index]; ++i_k_vector) {
                    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
                        for (size_t i_atom = 0; i_atom < atom_types_indexes[i_atom_type].size(); ++i_atom) {
                            unsigned int atom_index = atom_types_indexes[i_atom_type][i_atom];
                            double kr_vector_product = 0.0;
                            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                double delta_x = trajectory_[current_frame][atom_index][i_dimension] - trajectory_[initial_frame][atom_index][i_dimension];
                                kr_vector_product += k_vectors[k_index][i_k_vector][i_dimension] * delta_x;
                            }
                            Fs_kt_[k_index][time_point] += cos(k_min * kr_vector_product)*scattering_lengths_[i_atom_type]*scattering_lengths_[i_atom_type];
                        }
                    }
                }
            }
            Fs_kt_[k_index][time_point] *= normalization_factor[k_index];
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


void SelfIntermediateScattering::generate_k_vectors(unsigned int const & k_squared, vector< vector< unsigned int > > & k_vectors)
{
    if (dimension_ == 3) {
        //        if (method_of_k_generation_ == "moderate") {
        //            unsigned int mag_k_vector = static_cast< unsigned int > (sqrt(mag_kvec_sqr));
        //            for (unsigned int i = 0; i <= mag_k_vector; ++i) {
        //                unsigned int j_max = static_cast< unsigned int > (sqrt(mag_kvec_sqr - i*i));
        //                for (unsigned int j = 0; j <= j_max; ++j) {
        //                    unsigned int k_sqr = mag_kvec_sqr - i*i - j*j;
        //                    unsigned int k = static_cast< unsigned int > (sqrt(k_sqr));
        //                    // check if k is a perfect square integer
        //                    if ( fabs(k_sqr - k*k) < 1e-16) {
        //                        vector< unsigned int > inds = {i, j, k};
        //                        k_vectors.push_back( inds );
        //                    }
        //                }
        //            }
        //            return;
        //        }
        //        else if (method_of_k_generation_ == "fast") {
        unsigned int i_max = static_cast< unsigned int > (sqrt(k_squared/3.0));
        for (unsigned int i = 0; i <= i_max; ++i) {
            unsigned int j_max = static_cast< unsigned int > (sqrt((k_squared - i*i)/2.0));
            for (unsigned int j = i; j <= j_max; ++j) {
                unsigned int k_sqr = k_squared - i*i - j*j;
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
                        // Write down all remaining permuations
                        k_vectors.push_back( {i, k, j} );
                        k_vectors.push_back( {j, i, k} );
                        k_vectors.push_back( {j, k, i} );
                        k_vectors.push_back( {k, i, j} );
                        k_vectors.push_back( {k, j, i} );
                    }
                }
            }
        }
        //        }
    }
    else {
        vector< unsigned int > k_vector_temp(dimension_, 0);
        recurrsive_generate_k_vectors(k_squared, 0, k_vector_temp, k_vectors);
        
        return;
    }
    // add other sampling methods that are also generic in different dimensions
}


void SelfIntermediateScattering::recurrsive_generate_k_vectors(unsigned int const k_squared, unsigned int const dimension,
                                                    vector< unsigned int > & k_vector_temp, vector< vector< unsigned int > > & k_vectors)
{
    unsigned int sum_squares = 0;
    for (size_t i_dimension = 0; i_dimension < dimension; ++i_dimension) {
        sum_squares += k_vector_temp[i_dimension]*k_vector_temp[i_dimension];
    }
    unsigned int max_iterator = static_cast< unsigned int > (sqrt(k_squared - sum_squares));
    for (size_t i_k = 0; i_k <= max_iterator; ++i_k) {
        k_vector_temp[dimension] = i_k;
        if (dimension != dimension_ - 1) {
            recurrsive_generate_k_vectors(k_squared, dimension + 1, k_vector_temp, k_vectors);
        }
        else {
            unsigned int k_squared_current = 0;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                k_squared_current += k_vector_temp[i_dimension]*k_vector_temp[i_dimension];
            }
            if (k_squared_current == k_squared) {
                k_vectors.push_back(k_vector_temp);
            }
        }
    }
}


inline void SelfIntermediateScattering::compute_k_vectors(unsigned int & k_squared, vector< vector< unsigned int > > & k_vectors)
{
    generate_k_vectors( k_squared, k_vectors);
    unsigned int number_of_k_vectors = k_vectors.size();
    while (number_of_k_vectors == 0) {
        k_squared += 1;
        generate_k_vectors( k_squared, k_vectors);
        number_of_k_vectors = k_vectors.size();
    }
}


void SelfIntermediateScattering::write_Fs_kt()
{
    if (trajectory_delta_time_ == 0) {
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
	ofstream output_Fskt_file(output_file_name_);
	
	output_Fskt_file << "# Self intermediate scattering function for atom types: ";
    output_Fskt_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_Fskt_file << atom_types_[i_atom_type];
        output_Fskt_file << "(" << scattering_lengths_[i_atom_type] << ")";
        output_Fskt_file << " ";
    }
    output_Fskt_file << "}";
    output_Fskt_file << " in " << atom_group_ << endl;
    if (time_scale_type_ == "linear") {
        output_Fskt_file << "# using " << time_scale_type_ << " scale\n";
    }
    else if (time_scale_type_ == "log") {
        output_Fskt_file << "# using " << time_scale_type_ << " scale, logscale resulted in " << number_of_time_points_ - time_array_indexes_.size() << " repeated points ignored\n";
    }
    output_Fskt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    output_Fskt_file << "#" << endl;
    output_Fskt_file << "# k values" << endl;
    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
        output_Fskt_file << k_values_[k_index] << endl;
    }
    output_Fskt_file << "#" << endl;
    output_Fskt_file << "# t | Fs(k, t)" << endl;
    for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
        output_Fskt_file << time_array_indexes_[time_point]*trajectory_delta_time_;
        for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
            output_Fskt_file << "\t" << Fs_kt_[k_index][time_point];
        }
        output_Fskt_file << endl;
    }
    
    //    output_Fskt_file << "# ---------------------------" << endl;
    //    output_Fskt_file << "#               k_1st_row    " << endl;
    //    output_Fskt_file << "#            ----------------" << endl;
    //    output_Fskt_file << "# t_1st_col | Fs(k, t)_matrix" << endl;
    //    output_Fskt_file << "# ---------------------------" << endl;
    //
    //    output_Fskt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    //	output_Fskt_file << "            ";
    //    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
    //		output_Fskt_file << "\t" << k_values_[k_index];
    //	}
    //	for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
    //		output_Fskt_file << endl;
    //        output_Fskt_file << time_array_indexes_[time_point]*trajectory_delta_time_;
    //		for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
    //			output_Fskt_file << "\t" << Fs_kt_[k_index][time_point];
    //		}
    //	}
    
	output_Fskt_file.close();
}


// Function to check that all the parameters provided by the user for
// self intermediate scattering are usable.  Ensures that the code will
// not have with needing enough data points
void SelfIntermediateScattering::check_parameters() throw()
{
    if (number_of_frames_to_average_ > end_frame_ - start_frame_) {
        cerr << "ERROR: Cannot have the number of frames to average be greater than the number supplied\n";
        cerr << endl;
        exit(1);
    }
    
    if (end_frame_ == 0 && number_of_time_points_ == 0) {
        cerr << "ERROR: We require more information to proceed, either end_frame or number_of_time_points\n";
        cerr << "       must be povided for us to continue.";
        cerr << endl;
        exit(1);
    }
    
	if (time_scale_type_ == "linear") {
		if (end_frame_ == 0) {
			end_frame_ = start_frame_ + number_of_time_points_*frame_interval_ + number_of_frames_to_average_;
		}
		if (number_of_time_points_ == 0) {
            number_of_time_points_ = (end_frame_ - start_frame_ - number_of_frames_to_average_)/frame_interval_;
            
		}
		if (number_of_time_points_*frame_interval_ + number_of_frames_to_average_ > end_frame_ - start_frame_) {
			end_frame_ = start_frame_ + number_of_time_points_*frame_interval_ + number_of_frames_to_average_;
			cerr << "WARNING: the number of frames required is greater then the number supplied" << endl;
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
			end_frame_ = start_frame_ + pow(frame_interval_,number_of_time_points_)  + number_of_frames_to_average_;
		}
        if (number_of_time_points_ == 0) {
            number_of_time_points_ = static_cast<unsigned int>(log(end_frame_ - start_frame_ - number_of_frames_to_average_)/log(frame_interval_));
        }
        if (static_cast<unsigned int>(pow(frame_interval_, number_of_time_points_) + 0.5) + number_of_frames_to_average_ > end_frame_ - start_frame_) {
            end_frame_ = start_frame_ + static_cast<unsigned int>(pow(frame_interval_, number_of_time_points_) + 0.5)  + number_of_frames_to_average_;
			cerr << "WARNING: the number of frames required is greater then the number supplied\n";
			cerr << "   setting end frame to minimum value allowed: ";
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
	
	if (is_wrapped_) {
		cerr << "WARNING: the trajectory provided is not unwrapped\n";
		cerr << "       : We will unwrapp if for you, but user discretion\n";
		cerr << "       : is advised";
		cerr << endl;
	}
    
    if (atom_types_.empty()) {
        cerr << "ERROR: No atom types provided, nothing will be computed\n";
        cerr << endl;
        exit(1);
    }
    
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
    
    // check that output can be opened
    ofstream output_Fskt_file(output_file_name_);
    if (!output_Fskt_file) {
        cerr << "ERROR: Output file for self intermediate scattering function: "
        << "\033[1;25m"
        << output_file_name_
        << "\033[0m"
        << ", could not be opened"
        << endl;
        exit(1);
    }
    output_Fskt_file.close();
}


void SelfIntermediateScattering::compute_time_array()
{
    // first time point is always zero
    time_array_indexes_.push_back(0);
    
    double       total_time     = frame_interval_;
    unsigned int frame_previous = 0;
    
    // TODO switch to try catch
    for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {
        // check that we have not exceeded allowed frames
        assert(static_cast<unsigned int>(total_time) + number_of_frames_to_average_ < end_frame_ - start_frame_ && "Error: Not eneough frames for calculation on log time scale");
        
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


void SelfIntermediateScattering::determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                                        double & average_squared_scattering_length,
                                                        size_t & number_of_atoms)
{
    number_of_atoms = 0;
    
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        select_atoms(atom_types_indexes[i_atom_type], atom_types_[i_atom_type], atom_group_);
        
        number_of_atoms += atom_types_indexes[i_atom_type].size();
        // Generate average_scattering_length
        if (scattering_lengths_.size() != 0) {
            average_squared_scattering_length += atom_types_indexes[i_atom_type].size() * scattering_lengths_[i_atom_type] * scattering_lengths_[i_atom_type];
        }
    }
    
    // normalize the average scattering length
    if (scattering_lengths_.size() == 0) {
        average_squared_scattering_length = 1.0;
    }
    else {
        average_squared_scattering_length /= number_of_atoms;
    }
}


void SelfIntermediateScattering::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the self intermediate scattering function is: ";
    cout << status * 100.0/time_array_indexes_.size();
    cout << " \%";
    cout << flush;
}