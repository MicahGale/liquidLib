//
//  BondOrderParameter.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_BondOrderParameter_hpp
#define LiquidLib_BondOrderParameter_hpp

#include <vector>
#include <string>
#include <random>

#include "Trajectory.hpp"

using namespace std;

class BondOrderParameter : public Trajectory {
public:
    BondOrderParameter();
    virtual ~BondOrderParameter();
    
    void read_command_inputs(int argc, char * argv[]);
    void read_input_file();
    void compute_BOP();
    void write_BOP();
    
protected:
    
private:
    // private member functions
    void check_parameters() throw();
    void compute_harmonic(size_t const & atom1_index, size_t const & atom2_index, size_t const & current_frame,
                          size_t & number_neighbors, vector<double> & real_term, vector<double> & imaginary_term);
    void calculate_BOP(vector<double> const & real_term, vector<double> const & imaginary_term,
                       size_t const & number_neighbors, size_t const & time_point, size_t const & atom_index);
    void compute_time_array();
    void determine_atom_indexes(vector < vector < unsigned int > > & atom_types_indexes,
                                size_t & number_of_atoms);
    inline void print_status(size_t & status);
    
    // private member variables
    bool is_averaged_;
    bool add_bar_;
    
    string input_file_name_;
    string output_file_name_;
    string time_scale_type_;
    
    unsigned int number_of_time_points_;
    unsigned int bond_parameter_order_; // sets what order the bond order parameter uses, i.e. q_6, q_4....
    
    double frame_interval_;
    double max_cutoff_length_; // cut-off for determining the nearest neighbors
    
    vector< string > atom_types_;
    string atom_group_;
    
    vector< unsigned int > time_array_indexes_;
    vector< vector< double > > bond_order_parameter_; // first dimension is time point, second dimension is atom, or 0 if averaged
};

#endif /* defined(LiquidLib_BondOrderParameter_hpp) */
