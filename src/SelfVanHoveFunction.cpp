//
//  SelfVanHoveFunction.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "SelfVanHoveFunction.hpp"

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
#include <algorithm>

#include "Trajectory.hpp"

using namespace std;

SelfVanHoveFunction::SelfVanHoveFunction() :
	input_file_name_("Gs_rt.in"),
	output_file_name_("Gs_rt.txt"),
	atom_type_("all"),
	atom_group_("system"),
	time_scale_type_("linear"),
	number_of_bins_(200),
	number_of_time_points_(0),
	number_of_frames_to_average_(1),
	max_cutoff_length_(0.0),
	frame_interval_(1.0)
{
}


SelfVanHoveFunction::~SelfVanHoveFunction()
{
}


void SelfVanHoveFunction::read_command_inputs(int argc, char * argv[])
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


void SelfVanHoveFunction::read_input_file()
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
            if (output_file_name_ != "Gs_rt.txt") {
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
		if (input_word == "time_scale_tyoe") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			time_scale_type_ = input_word;
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
		if (input_word == "max_cutoff_length") {
			input_file >> input_word;
			if (input_word[0] == '=') {
				input_file >> input_word;
			}
			max_cutoff_length_ = stod(input_word);
			continue;
		}
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


void SelfVanHoveFunction::compute_Gs_rt()
{
	if (is_wrapped_) {
		unwrap_coordinates();
	}
    
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
	// select the indexes of atom_type_ in atom_group_
	vector< unsigned int > atom_type_indexes;
	select_atoms(atom_type_indexes, atom_type_, atom_group_);

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
    
    Gs_rt_.resize(number_of_bins_, vector< double >(time_array_indexes_.size(), 0.0));		// add one row Gs(r, t = 0)
    
    int status = 0;
    
	// Perform time averaging for Gs_rt_
#pragma omp parallel for
	for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {	// shift 1 to retain a row for Gs(r, t = 0)
		for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
			size_t current_frame = initial_frame + time_array_indexes_[time_point];
			for (vector< unsigned int >::iterator i_atom = atom_type_indexes.begin(); i_atom != atom_type_indexes.end(); ++i_atom) {
				double total_distance = 0.0;
				for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
					double delta_x = trajectory_[current_frame][*i_atom][i_dimension] - trajectory_[initial_frame][*i_atom][i_dimension];
					total_distance += delta_x * delta_x;
				}
				total_distance = sqrt(total_distance);
				
				unsigned int bin = round(total_distance/delta_r);
				if (bin < number_of_bins_) {
					Gs_rt_[bin][time_point] += 1.0;
				}
			}
		}
#pragma omp critical
{
        ++status;
        cout << "\rcurrent progress of calculating self van hove function is: ";
        cout << status * 100.0/number_of_time_points_;
        cout << " \%";
        cout << flush;
}
	}
	
	// do normalization
	r_values_.resize(number_of_bins_, 0.0);
	double dimension_scaling_factor = pow(M_PI, dimension_/2.0) / tgamma(1 + dimension_/2.0);	//Gamma function from <cmath>
	
	for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
		r_values_[i_bin] = delta_r * i_bin;
		
		// computer shell volume in a general dimension
		double volume_of_outer_sphere = pow(r_values_[i_bin] + delta_r/2, dimension_) * dimension_scaling_factor;
		double volume_of_inner_sphere = 0.0;
		if (i_bin != 0) {
		    volume_of_inner_sphere = pow(r_values_[i_bin] - delta_r/2, dimension_)* dimension_scaling_factor;
		}
		double volume_of_shell = volume_of_outer_sphere - volume_of_inner_sphere;
		
		double normalization_factor = 1.0/(volume_of_shell * atom_type_indexes.size() * number_of_frames_to_average_);
		
		if (i_bin == 0) {
			Gs_rt_[0][0] = 1.0/volume_of_shell;	// Gs(r, t = 0) = 1.0/volume_of_shell * delta_function(r)
		}
		for (size_t time_point = 1; time_point < number_of_time_points_; ++time_point) {
			Gs_rt_[i_bin][time_point] *= normalization_factor;
		}
	}
	
}


void SelfVanHoveFunction::write_Gs_rt()
{
    if (trajectory_delta_time_ == 0) {
        cerr << "WARNING: time step of simulation could not be derived from trajectory,\n";
        cerr << "       : and was not provided by input, will use time step of: ";
        cerr << "\033[1;25m" << "1 (step/a.u.)" << "\033[0m\n\n";
        trajectory_delta_time_ = 1.0;
    }
    
	ofstream output_Gsrt_file(output_file_name_);
	
	if (!output_Gsrt_file) {
		cerr << "ERROR: Output file for self van Hove function: "
			 << "\033[1;25m"
			 << output_file_name_
			 << "\033[0m"
			 << ", could not be opened"
			 << endl;
		exit(1);
	}
	
	output_Gsrt_file << "# Self van Hove function for atom type: ";
	output_Gsrt_file << atom_type_ << " in " << atom_group_ << endl;
	output_Gsrt_file << "# Data structure:" << endl;
	output_Gsrt_file << "# ---------------------------" << endl;
	output_Gsrt_file << "#               t_1st_row    " << endl;
	output_Gsrt_file << "#            ----------------" << endl;
	output_Gsrt_file << "# r_1st_col | Gs(r, t)_matrix" << endl;
	output_Gsrt_file << "# ---------------------------" << endl;
	
	output_Gsrt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
	output_Gsrt_file << "            ";
	for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
		output_Gsrt_file << "\t" << time_array_indexes_[time_point]*trajectory_delta_time_;
	}
	for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
		output_Gsrt_file << endl;
		output_Gsrt_file << r_values_[i_bin];
		for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
			output_Gsrt_file << "\t" << Gs_rt_[i_bin][time_point];
		}
	}
	
	output_Gsrt_file.close();
}


// Function to check that all the parameters provided by the user for
// self van Hove function are usable.  Ensures that the code will
// not have with needing enough data points
void SelfVanHoveFunction::check_parameters() throw()
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
}


void SelfVanHoveFunction::compute_time_array()
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
