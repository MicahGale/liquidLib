//
//  VelocityAutoCorrelationClass.hpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_VelocityAutoCorrelation_hpp
#define LiquidLib_VelocityAutoCorrelation_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class VelocityAutoCorrelation : public Trajectory {
public:
    VelocityAutoCorrelation();
    virtual ~VelocityAutoCorrelation();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_vacf_t();
    void write_vacf_t();

private:
//private member functions
    void check_parameters() throw();
    void compute_time_array();

// private variables
    string input_file_name_;
    string output_file_name_;
    string atom_type_;
    string atom_group_;
    string time_scale_type_;
    
    unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
    
    double frame_interval_;
    
    vector< unsigned int > time_array_indexes_;
    vector< double > vacf_t_;
};

#endif // defined (LiquidLib_VelocityAutoCorrelation_hpp)
