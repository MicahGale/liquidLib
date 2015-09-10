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
#ifndef LiquidLib_SelfIntermediateScattering_hpp
#define LiquidLib_SelfIntermediateScattering_hpp

#include <vector>
#include <string>
#include <random>

#include "Trajectory.hpp"

using namespace std;

class SelfIntermediateScattering : public Trajectory {
public:
    SelfIntermediateScattering();
    virtual ~SelfIntermediateScattering();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_Fs_kt();
    void write_Fs_kt();
   
protected:
    
private:
// private member functions
    void check_parameters() throw();
    void compute_time_array();
    void inline generate_k_vector(double const & k_absolute_value, vector< double > & k_vector);
    
// private member variables
    string input_file_name_;
    string output_file_name_;
    string atom_type_;
    string atom_group_;
	string time_scale_type_;
    string method_of_k_sampling_;       // analytical, gaussian, ...

    unsigned int number_of_bins_;       // bin spacing = 2*M_PI/average_box_length
    unsigned int k_start_index_;        // in unit of 2*M_PI/average_box_length
    unsigned int number_of_k_vectors_;
	unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
	
	double frame_interval_;
	
    vector< double > k_values_;
	vector< unsigned int > time_array_indexes_;
    vector< vector< double > > Fs_kt_;
    
    default_random_engine generator;
};

#endif // defined (LiquidLib_SelfIntermediateScattering_hpp)
