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
#include <sstream>

#include "Trajectory.hpp"

using namespace std;

SelfVanHoveFunction::SelfVanHoveFunction() :
	input_file_name_("Gs_rt.in"),
	output_file_name_("Gs_rt.txt"),
    atom_types_({}),
    scattering_lengths_({}),
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
			cerr << "ERROR: gro files cannot be used in non gromacs";
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


void SelfVanHoveFunction::compute_Gs_rt()
{
	if (is_wrapped_) {
		unwrap_coordinates();
	}
    
    // Form a array of time index values for a given type of timescale computation
    compute_time_array();
    
    // select the indexes of atom_types_
    size_t number_of_atoms;
    double average_squared_scattering_length = 0.0;
    vector < vector< unsigned int > > atom_types_indexes(atom_types_.size());
    
    determine_atom_indexes(atom_types_indexes, average_squared_scattering_length, number_of_atoms);

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
    
    size_t status = 1;
    cout << "Computing ..." << endl;
    
	// Perform time averaging for Gs_rt_
#pragma omp parallel for
	for (size_t time_point = 1; time_point < time_array_indexes_.size(); ++time_point) {	// shift 1 to retain a row for Gs(r, t = 0)
		for (size_t initial_frame = 0; initial_frame < number_of_frames_to_average_; ++initial_frame) {
			size_t current_frame = initial_frame + time_array_indexes_[time_point];
            for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
                for (size_t i_atom = 0; i_atom < atom_types_indexes.size(); ++i_atom) {
                    size_t atom_index = atom_types_indexes[i_atom_type][i_atom];
                    double total_distance = 0.0;
                    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                        double delta_x = trajectory_[current_frame][atom_index][i_dimension] - trajectory_[initial_frame][atom_index][i_dimension];
                        total_distance += delta_x * delta_x;
                    }
                    total_distance = sqrt(total_distance);
                    
                    unsigned int bin = round(total_distance/delta_r);
                    if (bin < number_of_bins_) {
                        Gs_rt_[bin][time_point] += 1.0;
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
    cout << endl;
	
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
		
		double normalization_factor = 1.0/(volume_of_shell * number_of_atoms * number_of_frames_to_average_);
        normalization_factor /= ( average_squared_scattering_length );
		
		if (i_bin == 0) {
			Gs_rt_[0][0] = 1.0/volume_of_shell;	// Gs(r, t = 0) = 1.0/volume_of_shell * delta_function(r)
		}
		for (size_t time_point = 1; time_point < time_array_indexes_.size(); ++time_point) {
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
	
	output_Gsrt_file << "# Self van Hove function for atom type: ";
    output_Gsrt_file << "# { ";
    for (size_t i_atom_type = 0; i_atom_type < atom_types_.size(); ++i_atom_type) {
        output_Gsrt_file << atom_types_[i_atom_type];
        output_Gsrt_file << "(" << scattering_lengths_[i_atom_type] << ")";
        output_Gsrt_file << " ";
    }
    output_Gsrt_file << "}";
    output_Gsrt_file << "in " << atom_group_ << ".\n";
    if (time_scale_type_ == "linear") {
        output_Gsrt_file << "# using " << time_scale_type_ << " scale\n";
    }
    else if (time_scale_type_ == "log") {
        output_Gsrt_file << "# using " << time_scale_type_ << " scale, logscale resulted in " << number_of_time_points_ - time_array_indexes_.size() << " repeated points ignored\n";
    }
    
    output_Gsrt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
    output_Gsrt_file << "#" << endl;
    output_Gsrt_file << "# t values" << endl;
    for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
        output_Gsrt_file << time_array_indexes_[time_point]*trajectory_delta_time_ << endl;
    }
    output_Gsrt_file << "#" << endl;
    output_Gsrt_file << "# r | Gs(r, t)" << endl;
    
    for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
        output_Gsrt_file << r_values_[i_bin];
        for (size_t time_point = 0; time_point < time_array_indexes_.size(); ++time_point) {
            output_Gsrt_file << "\t" << Gs_rt_[i_bin][time_point];
        }
        output_Gsrt_file << endl;
    }
    
//	output_Gsrt_file << "# ---------------------------" << endl;
//	output_Gsrt_file << "#               t_1st_row    " << endl;
//	output_Gsrt_file << "#            ----------------" << endl;
//	output_Gsrt_file << "# r_1st_col | Gs(r, t)_matrix" << endl;
//	output_Gsrt_file << "# ---------------------------" << endl;
//	
//	output_Gsrt_file << setiosflags(ios::scientific) << setprecision(output_precision_);
//	output_Gsrt_file << "            ";
//	for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
//		output_Gsrt_file << "\t" << time_array_indexes_[time_point]*trajectory_delta_time_;
//	}
//	for (size_t i_bin = 0; i_bin < number_of_bins_; ++i_bin) {
//		output_Gsrt_file << endl;
//		output_Gsrt_file << r_values_[i_bin];
//		for (size_t time_point = 0; time_point < number_of_time_points_; ++time_point) {
//			output_Gsrt_file << "\t" << Gs_rt_[i_bin][time_point];
//		}
//	}
	
	output_Gsrt_file.close();
}


// Function to check that all the parameters provided by the user for
// self van Hove function are usable.  Ensures that the code will
// not have with needing enough data points
void SelfVanHoveFunction::check_parameters() throw()
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
    output_Gsrt_file.close();
}


void SelfVanHoveFunction::compute_time_array()
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


void SelfVanHoveFunction::determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
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


void SelfVanHoveFunction::print_status(size_t & status)
{
    ++status;
    cout << "\rcurrent progress of calculating the self van hove function is: ";
    cout << status * 100.0/time_array_indexes_.size();
    cout << " \%";
    cout << flush;
}