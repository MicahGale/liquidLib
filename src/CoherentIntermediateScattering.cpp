//
//  CoherentIntermediateScattering.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "CoherentIntermediateScattering.hpp"

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

#include "Trajectory.hpp"

using namespace std;

CoherentIntermediateScattering::CoherentIntermediateScattering() :
	input_file_name_("F_kt.in"),
	output_file_name_("F_kt.txt"),
	atom_group_("system"),
	time_scale_type_("linear"),
    method_of_k_generation_("moderate"),
    number_of_bins_(50),
    k_start_index_(1),
	number_of_time_points_(0),
	number_of_frames_to_average_(1),
	frame_interval_(1.0),
    delta_k_(0)
{
}


CoherentIntermediateScattering::~CoherentIntermediateScattering()
{
}


void CoherentIntermediateScattering::read_command_inputs(int argc, char * argv[])
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


void CoherentIntermediateScattering::read_input_file()
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
            if (output_file_name_ != "F_kt.txt") {
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
        if (input_word == "method_of_k_generation") {
            input_file >> input_word;
            if (input_word[0] == '=') {
                input_file >> input_word;
            }
            method_of_k_generation_ = input_word;
            cout << method_of_k_generation_ << endl;
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


void CoherentIntermediateScattering::compute_F_kt()
{
	if (is_wrapped_) {
		unwrap_coordinates();
	}
	
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
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
    double k_min = 2.0*M_PI/min_box_length;
    delta_k_ = (delta_k_ < k_min ) ? k_min : delta_k_;
    
    // allocate k_values_ and normalization factor
    k_values_.resize(number_of_bins_, 0.0);
    number_of_k_vectors_.resize(number_of_bins_);
    vector< vector< vector < unsigned int > > > k_vectors(number_of_bins_);
    vector< double > normalization_factor(number_of_bins_, 0.0);
    
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        
        // Compare delta_k_ with k_min allowed to select the larger of the two for first k-point
        double k_value_temp = delta_k_ * (k_index + k_start_index_);
        
        // Compute the integer value based on user input k_value
        unsigned int mag_kvec_sqr = static_cast<unsigned int> (k_value_temp * k_value_temp * min_box_length * min_box_length / (4.0 * M_PI * M_PI));
        
        // Compute the avaible k_vectors
        compute_k_vectors( mag_kvec_sqr, k_vectors[k_index] );
        number_of_k_vectors_[k_index] = k_vectors[k_index].size();
        
        // Corrected k_value on the grid with the new mag_kvec_sqr
        k_values_[k_index] = 2.0 * M_PI * sqrt(1.0 * mag_kvec_sqr) / min_box_length;
        
        // Compute normalization factor
        normalization_factor[k_index] = 1.0 / (number_of_frames_to_average_ * number_of_atoms * number_of_k_vectors_[k_index]);
        normalization_factor[k_index] /= ( average_scattering_length * average_scattering_length );
        
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
    cout << "Finished computing correct k-values and vectors. " << endl;
    
    // initialize F_kt_
    F_kt_.resize(number_of_bins_, vector< double > (time_array_indexes_.size(), 0.0));

    // Store periodic_image_counts of trajectories until number_of_frames_to_average
    double *** periodic_image_counts;
    periodic_image_counts = new double ** [number_of_frames_to_average_];
#pragma omp parallel for
    for (size_t i_frame = 0; i_frame < number_of_frames_to_average_; ++i_frame) {
        periodic_image_counts[i_frame] = new double * [number_of_system_atoms_];
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            periodic_image_counts[i_frame][i_atom] = new double [dimension_];
        }
    }
    for (size_t i_frame = 0; i_frame < number_of_frames_to_average_; ++i_frame) {
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                periodic_image_counts[i_frame][i_atom][i_dimension] = box_length_[i_frame][i_dimension] * round( trajectory_[i_frame][i_atom][i_dimension] / box_length_[i_frame][i_dimension] );
            }
        }
    }
    
    size_t status = 0;
    cout << "Computing ..." << endl;
    
#pragma omp parallel for firstprivate(k_vector, k_vectors)
    for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
        for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
            
            for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
                size_t current_frame = initial_frame + time_array_indexes_[time_point];
                
                for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_[k_index]; ++i_k_vector) {
                    double sum_cos_term_r0 = 0.0;
                    double sum_sin_term_c0 = 0.0;
                    double sum_cos_term_rt = 0.0;
                    double sum_sin_term_ct = 0.0;
                    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
                        for (size_t i_atom = 0; i_atom < atom_types_indexes[i_atom_type].size(); ++i_atom) {
                            size_t atom_index = atom_types_indexes[i_atom_type][i_atom];
                            double k_dot_r0 = 0.0;
                            double k_dot_rt = 0.0;
                            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                k_dot_r0 += k_vectors[k_index][i_k_vector][i_dimension] * (trajectory_[initial_frame][atom_index][i_dimension] - periodic_image_counts[initial_frame][atom_index][i_dimension]);
                                k_dot_rt += k_vectors[k_index][i_k_vector][i_dimension] * (trajectory_[current_frame][atom_index][i_dimension] - periodic_image_counts[initial_frame][atom_index][i_dimension]);
                            }
                            sum_cos_term_r0 += cos(k_min * k_dot_r0);
                            sum_sin_term_c0 += sin(k_min * k_dot_r0);
                            sum_cos_term_rt += cos(k_min * k_dot_rt);
                            sum_sin_term_ct += sin(k_min * k_dot_rt);
                        }
                        sum_cos_term_r0 *= scattering_lengths_[i_atom_type];
                        sum_sin_term_c0 *= scattering_lengths_[i_atom_type];
                        sum_cos_term_rt *= scattering_lengths_[i_atom_type];
                        sum_sin_term_ct *= scattering_lengths_[i_atom_type];
                    }
                    F_kt_[k_index][time_point] += (sum_cos_term_rt * sum_cos_term_r0 + sum_sin_term_ct * sum_sin_term_c0);
                }
            }
            F_kt_[k_index][time_point] *= normalization_factor[k_index];
        }
        
        if (is_run_mode_verbose_) {
#pragma omp atomic
            print_status(status);
        }
    }
    cout << endl;
    
    // delete dynamically allocated memory
    for (size_t i_frame = 0; i_frame < number_of_frames_to_average_; ++i_frame) {
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            delete[] periodic_image_counts[i_frame][i_atom];
        }
        delete[] periodic_image_counts[i_frame];
    }
    delete[] periodic_image_counts;
}


void CoherentIntermediateScattering::generate_k_vectors(unsigned int const & mag_kvec_sqr, vector< vector< unsigned int > > & k_vectors)
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
                    if (i != j && i != k && j != k) {
                        // Write down all 6 permuations
                        k_vectors.push_back( {i, k, j} );
                        k_vectors.push_back( {j, i, k} );
                        k_vectors.push_back( {j, k, i} );
                        k_vectors.push_back( {k, i, j} );
                        k_vectors.push_back( {k, j, i} );
                    }
                    else if ( i == j || i == k || j == k ) {
                        k_vectors.push_back( {j, k, i} );
                        k_vectors.push_back( {k, i, j} );
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


void CoherentIntermediateScattering::recurrsive_k_generation(unsigned int const mag_kvec_sqr, unsigned int const dimension,
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


inline void CoherentIntermediateScattering::compute_k_vectors(unsigned int & mag_kvec_sqr, vector< vector< unsigned int > > & k_vectors)
{
    generate_k_vectors( mag_kvec_sqr, k_vectors);
    unsigned int number_of_k_vectors = k_vectors.size();
    while (number_of_k_vectors == 0) {
        mag_kvec_sqr += 1;
        generate_k_vectors( mag_kvec_sqr, k_vectors);
        number_of_k_vectors = k_vectors.size();
    }
}


void CoherentIntermediateScattering::write_F_kt()
{
    if (trajectory_delta_time_ == 0) {
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
	ofstream output_Fkt_file(output_file_name_);
	
	output_Fkt_file << "# Coherent intermediate scattering function for atom type: \n";
    output_Fkt_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_Fkt_file << atom_types_[i_atom_type];
        output_Fkt_file << "(" << scattering_lengths_[i_atom_type] << ")";
        output_Fkt_file << "  ";
    }
    output_Fkt_file << "}";
	output_Fkt_file << " in " << atom_group_ << endl;
    output_Fkt_file << "# using " << time_scale_type_ << " time scale, " << endl;
    
    output_Fkt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    output_Fkt_file << "#" << endl;
    output_Fkt_file << "# k values" << endl;
    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
        output_Fkt_file << k_values_[k_index] << endl;
    }
    output_Fkt_file << "#" << endl;
    output_Fkt_file << "# t | F(k, t)" << endl;
    for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
        output_Fkt_file << time_array_indexes_[time_point]*trajectory_delta_time_;
        for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
            output_Fkt_file << "\t" << F_kt_[k_index][time_point];
        }
        output_Fkt_file << endl;
    }
    
    output_Fkt_file.close();
}


// Function to check that all the parameters provided by the user for
// Coherent intermediate scattering are usable.  Ensures that the code will
// not have with needing enough data points
void CoherentIntermediateScattering::check_parameters() throw()
{
    if (end_frame_ == 0 && number_of_time_points_ == 0) {
        cerr << "ERROR: We require more information to proceed, either frameend or numberoftimepoints\n";
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
	
	if (is_wrapped_) {
		cerr << "WARNING: the trajectory provided is not unwrapped\n";
		cerr << "       : We will unwrapp if for you, but user discretion\n";
		cerr << "       : is advised";
		cerr << endl;
	}
    
    if (k_start_index_ < 1) {
        cerr << "WARNING: k_start_index cannot be non-positive\n";
        cerr << "       : We will reset the value to 1\n";
        cerr << endl;
        k_start_index_ = 1;
    }
    
    if (atom_types_.empty()) {
        cerr << "\n";
        cerr << "ERROR: No atom types specified.\n";
        cerr << "     : Computation cannot proceed.\n";
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
    
    // check output file can be opened
    ofstream output_Fkt_file(output_file_name_);
    if (!output_Fkt_file) {
        cerr << "ERROR: Output file for Coherent intermediate scattering function: "
        << "\033[1;25m"
        << output_file_name_
        << "\033[0m"
        << ", could not be opened"
        << endl;
        exit(1);
    }
    output_Fkt_file.close();
}


void CoherentIntermediateScattering::determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
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


void CoherentIntermediateScattering::compute_time_array()
{
    time_array_indexes_.resize(number_of_time_points_);
    time_array_indexes_[0] = 0;
    
    double       total_time     = frame_interval_;
    unsigned int frame_previous = 0;
    
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
        
        assert(time_array_indexes_[time_point] + number_of_frames_to_average_ < end_frame_ - start_frame_ && "Error: Not eneough frames for calculation on log time scale");
    }
}


void CoherentIntermediateScattering::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the pair distribution function is: ";
    cout << status * 100.0/number_of_time_points_;
    cout << " \%";
    cout << flush;
}