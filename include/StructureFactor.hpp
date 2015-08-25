//
//  StructureFactor.hpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_StructureFactor_hpp
#define LiquidLib_StructureFactor_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class StructureFactor : public Trajectory {
public:
    StructureFactor();
    virtual ~StructureFactor();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_S_k();
    void write_S_k();
    
protected:
    
private:
// private member functions
    void check_parameters() throw();
    void generate_k_vectors(vector< vector< double > > & k_vectors, double const & k_absolute_value);
	
// private member variables
    string input_file_name_;
    string output_file_name_;
    string atom_group_;
    string method_of_k_sampling_;       // gaussian, ...

    unsigned int number_of_bins_;       // bin spacing = 2*M_PI/average_box_length
    unsigned int k_start_index_;        // in unit of 2*M_PI/average_box_length
    unsigned int number_of_k_vectors_;  // for k sampling, ingoring analytical method
    unsigned int number_of_frames_to_average_;
	
    vector< string > user_atom_types_;
    
    vector< double > k_values_;
    vector< double > scattering_lengths_;
    vector< double > S_k_;
};

#endif // defined (LiquidLib_StructureFactor_hpp)
