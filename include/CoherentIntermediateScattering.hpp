//
//  CoherentIntermediateScattering.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_CoherentIntermediateScattering_hpp
#define LiquidLib_CoherentIntermediateScattering_hpp

#include <vector>
#include <string>
#include <random>

#include "Trajectory.hpp"

using namespace std;

class CoherentIntermediateScattering : public Trajectory {
public:
    CoherentIntermediateScattering();
    virtual ~CoherentIntermediateScattering();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_F_kt();
    void write_F_kt();
   
protected:
    
private:
// private member functions
    void check_parameters() throw();
    void compute_time_array();
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                double & average_scattering_length,
                                size_t & number_of_atoms);
    void generate_k_vectors(unsigned int const & mag_kvec_sqr, vector< vector< unsigned int > > & kvec_temp);
    inline void compute_k_vectors(unsigned int & mag_kvec_sqr, vector< vector< unsigned int > > & kvec_temp);
    inline void print_status(size_t & status);
    
// private member variables
    string input_file_name_;
    string output_file_name_;
    string atom_group_;
	string time_scale_type_;
    string method_of_k_generation_;

    unsigned int number_of_bins_;       // bin spacing = 2*M_PI/average_box_length
    unsigned int k_start_index_;
	unsigned int number_of_time_points_;
    unsigned int number_of_frames_to_average_;
	
	double frame_interval_;
    double delta_k_;
    
    vector< string > atom_types_;
	
    vector< double > k_values_;
    vector< double > scattering_lengths_;
	vector< unsigned int > time_array_indexes_;
    vector< unsigned int > number_of_k_vectors_;
    vector< vector< double > > F_kt_;
};

#endif // defined (LiquidLib_CoherentIntermediateScattering_hpp)
