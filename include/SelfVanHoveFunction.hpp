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
#ifndef LiquidLib_SelfVanHoveFunction_hpp
#define LiquidLib_SelfVanHoveFunction_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class SelfVanHoveFunction : public Trajectory {
public:
    SelfVanHoveFunction();
    virtual ~SelfVanHoveFunction();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_Gs_rt();
    void write_Gs_rt();
  
protected:
    
private:
// private member functions
    void check_parameters() throw();
    void compute_time_array();
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                double & average_scattering_length,
                                size_t & number_of_atoms);
    inline void print_status(size_t & status);
	
// private member variables
    string input_file_name_;
    string output_file_name_;
    vector < string > atom_types_;
    vector < double > scattering_lengths_;
    string atom_group_;
	string time_scale_type_;

    unsigned int number_of_bins_;
	unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
	
	double max_cutoff_length_;
	double frame_interval_;
	
    vector< double > r_values_;
	vector< unsigned int > time_array_indexes_;
    vector< vector< double > > Gs_rt_;
};

#endif // defined (LiquidLib_SelfVanHoveFunction_hpp)
