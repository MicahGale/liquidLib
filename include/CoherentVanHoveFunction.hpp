//
//  CoherentVanHoveFunction.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_CoherentVanHoveFunction_hpp
#define LiquidLib_CoherentVanHoveFunction_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class CoherentVanHoveFunction : public Trajectory {
public:
    CoherentVanHoveFunction();
    virtual ~CoherentVanHoveFunction();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_G_rt();
    void write_G_rt();
  
protected:
    
private:
// private member functions
    void check_parameters() throw();
    void compute_time_array();
	
// private member variables
    string input_file_name_;
    string output_file_name_;
	string time_scale_type_;
    string atom_group_;
    vector< string > user_atom_types_;
    vector< double > scattering_lengths_;
    
    unsigned int number_of_bins_;
	unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
	
	double max_cutoff_length_;
	double frame_interval_;
	
    vector< double > r_values_;
	vector< unsigned int > time_array_indexes_;
    vector< vector< double > > G_rt_;
};

#endif // defined (LiquidLib_CoherentVanHoveFunction_hpp)
