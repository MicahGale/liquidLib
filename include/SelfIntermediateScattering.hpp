//
//  SelfIntermediateScattering.cpp
//
//  Copyright 2015 Zhang-Group, 2020 Micah Gale.
//  All rights reserved.
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
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                double & average_squared_scattering_length,
                                size_t & number_of_atoms);
    void generate_k_vectors(unsigned int const & k_squared, vector< vector< unsigned int > > & k_vectors);
    void recurrsive_generate_k_vectors(unsigned int const k_squared, unsigned int const dimension,
                                       vector< unsigned int > & k_vector_temp, vector< vector< unsigned int > > & k_vectors);
    inline void compute_k_vectors(unsigned int & k_squared, vector< vector< unsigned int > > & k_vectors);
    inline void print_status(size_t & status);
    
// private member variables
    string input_file_name_;
    string output_file_name_;
    string atom_group_;
	string time_scale_type_;

    unsigned int number_of_bins_;       // bin spacing = 2*M_PI/average_box_length
    unsigned int k_start_index_;        // in unit of 2*M_PI/average_box_length
	unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
    unsigned int number_of_k_;          // the number of input k-values
	
	double frame_interval_;
    double delta_k_;
    
    vector< string > atom_types_;
	
    vector< double > k_values_;
    vector< double > scattering_lengths_;
    vector< double > input_k_values_;
	vector< unsigned int > time_array_indexes_;
    vector< unsigned int > number_of_k_vectors_;
    vector< vector< double > > Fs_kt_;
};

#endif // defined (LiquidLib_SelfIntermediateScattering_hpp)
