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

#include "Trajectory.hpp"

using namespace std;

SelfIntermediateScattering::SelfIntermediateScattering() :
	input_file_name_("Fs_kt.in"),
	output_file_name_("Fs_kt.txt"),
	atom_type_("all"),
	atom_group_("system"),
	time_scale_type_("linear"),
    method_of_k_sampling_("analytical"),
    number_of_bins_(50),
    k_start_index_(0),
    number_of_k_vectors_(50),
	number_of_time_points_(0),
	number_of_frames_to_average_(1),
	frame_interval_(1.0)
{
}


SelfIntermediateScattering::~SelfIntermediateScattering()
{
}


void SelfIntermediateScattering::read_command_inputs(int argc, char * argv[])
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


void SelfIntermediateScattering::compute_Fs_kt()
{
	if (is_wrapped_) {
		unwrap_coordinates();
	}
	
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
	// select the indexes of atom_type_ in atom_group_
	vector< unsigned int > atom_type_indexes;
	select_atoms(atom_type_indexes, atom_type_, atom_group_);

    // k resolution is determined by inverse box length
    double min_box_length = average_box_length_[0];
    for (size_t i_dimension = 1; i_dimension < dimension_; ++i_dimension) {
        min_box_length = (min_box_length < average_box_length_[i_dimension]) ? min_box_length : average_box_length_[i_dimension];
    }
    double delta_k = 2.0*M_PI/min_box_length;

    // allocate k_values_ and Fs_kt_
    k_values_.resize(number_of_bins_, 0.0);
    for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
        k_values_[k_index] = delta_k * (k_index + k_start_index_);
    }
    Fs_kt_.resize(number_of_bins_, vector< double >(time_array_indexes_.size(), 1.0));   // initialize with value 1.0 since Fs(k, t = 0) = 1.0

    double normalization_factor = 1.0/(number_of_frames_to_average_ * atom_type_indexes.size() * number_of_k_vectors_);
    
    int status = 0;
    
    vector< vector < double > > k_vectors(number_of_k_vectors_, vector< double >(dimension_, 0.0));
    
#pragma omp parallel for firstprivate(k_vectors)
    for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {
        for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
            double sum_of_individual_terms = 0.0;
            
            for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
                size_t current_frame = initial_frame + time_array_indexes_[time_point];

                if (method_of_k_sampling_ == "analytical") {
                    for (vector< unsigned int >::iterator i_atom = atom_type_indexes.begin(); i_atom != atom_type_indexes.end(); ++i_atom) {
                        double total_distance = 0.0;
                        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                            double delta_x = trajectory_[current_frame][*i_atom][i_dimension] - trajectory_[initial_frame][*i_atom][i_dimension];
                            total_distance += delta_x * delta_x;
                        }
                        total_distance = sqrt(total_distance);
                        double k_times_r = k_values_[k_index] * total_distance;
                        sum_of_individual_terms += sin(k_times_r)/k_times_r;
                    }
                }
                else {
                    generate_k_vectors(k_vectors, k_values_[k_index]);
                    for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_; ++i_k_vector) {
                        for (vector< unsigned int >::iterator i_atom = atom_type_indexes.begin(); i_atom != atom_type_indexes.end(); ++i_atom) {
                            double kr_vector_product = 0.0;
                            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                                double delta_x = trajectory_[current_frame][*i_atom][i_dimension] - trajectory_[initial_frame][*i_atom][i_dimension];
                                kr_vector_product += k_vectors[i_k_vector][i_dimension] * delta_x;
                            }
                            sum_of_individual_terms += cos(kr_vector_product);
                        }
                    }
                }
                
            }
            Fs_kt_[k_index][time_point] = sum_of_individual_terms * normalization_factor;
        }
#pragma omp critical
{
        ++status;
        cout << "\rcurrent progress of calculating self intermediate scattering is: ";
        cout << status * 100.0/number_of_time_points_;
        cout << " \%";
        cout << flush;
}
    }
	
}


void SelfIntermediateScattering::generate_k_vectors(vector< vector< double > > & k_vectors, double const & k_absolute_value)
{
    if (method_of_k_sampling_ == "gaussian") {
        random_device seed;
        default_random_engine generator(seed());
        normal_distribution< double > distribution(0.0, 1.0);
        vector< double > random_vector(dimension_, 0.0);
        for (size_t i_k_vector = 0; i_k_vector < number_of_k_vectors_; ++i_k_vector) {
            double vector_length = 0.0;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                random_vector[i_dimension] = distribution(generator);
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


void SelfIntermediateScattering::write_Fs_kt()
{
    if (trajectory_delta_time_ == 0) {
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
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
	
	output_Fskt_file << "# Self intermediate scattering function for atom type: ";
	output_Fskt_file << atom_type_ << " in " << atom_group_ << endl;
    output_Fskt_file << "# using " << time_scale_type_ << " time scale, ";
    output_Fskt_file << method_of_k_sampling_ << " k sampling" << endl;
    
//	output_Fskt_file << "# Data structure:" << endl;
//	output_Fskt_file << "# ---------------------------" << endl;
//	output_Fskt_file << "#               t_1st_row    " << endl;
//	output_Fskt_file << "#            ----------------" << endl;
//	output_Fskt_file << "# k_1st_col | Fs(k, t)_matrix" << endl;
//	output_Fskt_file << "# ---------------------------" << endl;
//	
//	output_Fskt_file << setiosflags(ios::scientific) << setprecision(6);
//	output_Fskt_file << "            ";
//	for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
//		output_Fskt_file << "\t" << time_array_indexes_[time_point]*trajectory_delta_time_;
//	}
//	for (size_t k_index = 0; k_index < number_of_bins_; ++k_index) {
//		output_Fskt_file << endl;
//		output_Fskt_file << k_values_[k_index];
//		for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
//			output_Fskt_file << "\t" << Fs_kt_[k_index][time_point];
//		}
//	}
    
    output_Fskt_file << "# Data structure:" << endl;
    output_Fskt_file << "# ---------------------------" << endl;
    output_Fskt_file << "#               k_1st_row    " << endl;
    output_Fskt_file << "#            ----------------" << endl;
    output_Fskt_file << "# t_1st_col | Fs(k, t)_matrix" << endl;
    output_Fskt_file << "# ---------------------------" << endl;
    
    output_Fskt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
	output_Fskt_file << "            ";
    for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
		output_Fskt_file << "\t" << k_values_[k_index];
	}
	for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
		output_Fskt_file << endl;
        output_Fskt_file << time_array_indexes_[time_point]*trajectory_delta_time_;
		for (size_t k_index = 0; k_index < k_values_.size(); ++k_index) {
			output_Fskt_file << "\t" << Fs_kt_[k_index][time_point];
		}
	}
    
	output_Fskt_file.close();
}


// Function to check that all the parameters provided by the user for
// self intermediate scattering are usable.  Ensures that the code will
// not have with needing enough data points
void SelfIntermediateScattering::check_parameters() throw()
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
		if (number_of_time_points_*frame_interval_ + number_of_frames_to_average_ > end_frame_ - start_frame_) {
			end_frame_ = start_frame_ + pow(frame_interval_,number_of_time_points_)  + number_of_frames_to_average_;
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
    
    if (method_of_k_sampling_ == "analytical" && number_of_k_vectors_ != 1) {
        number_of_k_vectors_ = 1;
        cout << "Note: Parameter \"number of k vectors\" is not needed for \"analytical\" k sampling and will be ignored" << endl;
        cout << endl;
    }
    
    if (method_of_k_sampling_ != "analytical" && method_of_k_sampling_ != "gaussian" && method_of_k_sampling_ != "uniform") {
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


void SelfIntermediateScattering::compute_time_array()
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
