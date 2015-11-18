//
//  MeanSquaredDisplacementClass.hpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_MeanSquaredDisplacement_hpp
#define LiquidLib_MeanSquaredDisplacement_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class MeanSquaredDisplacement : public Trajectory {
public:
    MeanSquaredDisplacement();
    virtual ~MeanSquaredDisplacement();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_r2_t();
    void write_r2_t();
    
protected:
//protected member functions
    void check_parameters() throw();
    void compute_time_array();
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                size_t & number_of_atoms);
    inline void print_status(size_t & status);
    
//protected variables
    unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
    
    double frame_interval_;
    
    vector< unsigned int > time_array_indexes_;
    vector< double > r2_t_;
    
    string input_file_name_;
    string output_file_name_;
    string time_scale_type_;
    vector < string > atom_types_;
    string atom_group_;
    
private:
};

#endif // defined (LiquidLib_MeanSquaredDisplacement_hpp)
