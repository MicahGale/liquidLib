//
//  Trajectory.hpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_Trajectory_hpp
#define LiquidLib_Trajectory_hpp

#ifdef GROMACS
extern "C" {
    #include "xdrfile.h"
}
#endif

#include <vector>
#include <string>

using namespace std;

class Trajectory {
public:
    Trajectory();
    virtual ~Trajectory();
    void read_trajectory();
    virtual void read_input_file() = 0;
    
protected:
// protected member functions
    void wrap_coordinates();
    void unwrap_coordinates();
    void compute_velocities();
    void select_atoms(vector< unsigned int > & selected_atom_indexes,
                      string const & atom_type_to_select = "all",
                      string const & atom_group_to_select = "system");
    
// protected variables
    bool is_run_mode_verbose_;                // sets if run output should be verbose
    bool is_wrapped_;
    
    unsigned int start_frame_;
    unsigned int end_frame_;
    unsigned int number_of_system_atoms_;
    unsigned int dimension_;
    
    vector< unsigned int > number_of_atoms_of_each_type_;
    
    double trajectory_delta_time_;
    double output_precision_;
    
    string trajectory_file_name_;
    string trajectory_file_type_;
    string gro_file_name_;
    string trajectory_data_type_;    // coordinate, velocity, or force

    vector< string > system_atom_types_; //contains atom_type of every atom
    vector< string > molecule_id_;
    
    // box_length_[time][coordinate]
    vector< double > average_box_length_;
    vector< vector< double > > box_length_;

    // trajectory_[time][atom][coordinate]
    double *** trajectory_;
    
private:
// private member functions
#ifdef GROMACS
    void read_gro_file();
	void read_xtc_file();
    void read_trr_file();
#endif
    
    void read_vasp_file();
    void read_dump_file();
    void read_abc_file();
    
#ifdef GROMACS
    inline void save_frame(unsigned int const & frame_number, float const & time, matrix const & box, rvec const * coordinate);
    inline void set_monatomic();
#endif
    
    inline void allocate_memory();
};

#endif // defined(LiquidLib_Trajectory_hpp)
