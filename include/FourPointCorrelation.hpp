//
//  FourPointCorrelation.hpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_FourPointCorrelation_hpp
#define LiquidLib_FourPointCorrelation_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class FourPointCorrelation : public Trajectory {
public:
    FourPointCorrelation();
    virtual ~FourPointCorrelation();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_chi4_t();
    void write_chi4_t();

protected:
    
private:
// private member functions
    void check_parameters() throw();
    void compute_time_array();
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                size_t & number_of_atoms);
    inline void print_status(size_t & status);
    
//private variables
    unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
    
    double frame_interval_;
    double overlap_length_;
    
    vector< vector< double > > chi4_t_;
    vector< unsigned int > time_array_indexes_;
    
    string time_scale_type_;
    string output_file_name_;
    string input_file_name_;
    vector < string > atom_types_;
    string atom_group_;
};

#endif // defined (LiquidLib_FourPointCorrelation_hpp)
