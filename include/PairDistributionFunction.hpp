//
//  PairDistributionFunction.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_PairDistributionFuntion_hpp
#define LiquidLib_PairDistributionFuntion_hpp

#include <vector>
#include <string>

#include "Trajectory.hpp"

using namespace std;

class PairDistributionFunction : public Trajectory {
public:
    PairDistributionFunction();
    virtual ~PairDistributionFunction();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_g_r();
    void write_g_r();

protected:
    
private:
// private member function
    void check_parameters() throw();
    inline void histogram_g_r(size_t const & frame_number, size_t const & i_atom1, size_t const & i_atom2, double const & delta_r);
    
// private member variables
    string input_file_name_;
    string output_file_name_;
    vector < string > atom_types1_;
    vector < string > atom_types2_;
    string atom_group1_;
    string atom_group2_;

    unsigned int number_of_bins_;
    unsigned int number_of_frames_to_average_;
	
	double max_cutoff_length_;
	
    vector< double > r_values_;
    vector< double > g_r_;
};

#endif // defined (LiquidLib_PairDistributionFuntion_hpp)
