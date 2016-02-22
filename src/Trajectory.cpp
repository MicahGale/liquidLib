//
//  Trajectory.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "Trajectory.hpp"

#ifdef OMP
#include "omp.h"
#endif

#ifdef GROMACS
extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_trr.h"
    #include "xdrfile_xtc.h"
}
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>

using namespace std;

// default constructor
Trajectory::Trajectory():
    is_run_mode_verbose_(false),
    is_wrapped_(false),
    start_frame_(0),
    end_frame_(0),
    number_of_system_atoms_(0),
    dimension_(3),
    trajectory_delta_time_(0.0),
    output_precision_(15),
    trajectory_file_name_(""),
    trajectory_file_type_(""),
    gro_file_name_(""),
    trajectory_data_type_("coordinate"),
    trajectory_(0)
{
}


// deconstructor
Trajectory::~Trajectory()
{
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            delete[] trajectory_[i_frame][i_atom];
        }
        delete[] trajectory_[i_frame];
    }
    delete[] trajectory_;
}


// function to determine file type of trajectory
// and then call corresponding read trajectory function
void Trajectory::read_trajectory()
{
    // Be sure that the trajectory file exists
    if (!ifstream(trajectory_file_name_)) {
        cerr << "ERROR: The trajectory file:";
        cerr << "\033[1;34m";
        cerr << trajectory_file_name_;
        cerr << "\033[0m";
        cerr << " does not exist.";
        cerr << endl;
        exit(1);
    }
#ifdef GROMACS
    // Be sure that the .gro file exists
    if (!ifstream(gro_file_name_) && gro_file_name_ != "") {
        cerr << "ERROR: The gro file:";
        cerr << "\033[1;34m";
        cerr << gro_file_name_;
        cerr << "\033[0m";
        cerr << " does not exists";
        cerr << endl;
        exit(1);
    }
#endif

#ifdef GROMACS
    if (trajectory_file_name_.rfind(".xtc") != string::npos) {
        read_xtc_file();
        if (gro_file_name_ != "") {
            read_gro_file();
        }
        else {
            set_monatomic();
        }        return;
    }
    if (trajectory_file_name_.rfind(".trr") != string::npos) {
        read_trr_file();
        if (gro_file_name_ != "") {
            read_gro_file();
        }
        else {
            set_monatomic();
        }
        return;
    }
#else
    if (trajectory_file_name_.rfind(".trr") != string::npos) {
        cerr << "ERROR: .trr files are not compatible in non GROMACS compatible version\n";
        cerr << "     : of LiquidLib.  Please add appropriate flags in Makefile\n";
        cerr << "     : and recompile the binaries.\n";
        cerr << endl;
        exit(1);
    }
    if (trajectory_file_name_.rfind(".xtc") != string::npos) {
        cerr << "ERROR: .xtc files are not compatible in non GROMACS version\n";
        cerr << "     : of LiquidLib.  Please add appropriate flags in Makefile\n";
        cerr << "     : and recompile the binaries.\n";
        cerr << endl;
        exit(1);
    }
#endif
    
    if (trajectory_file_name_.rfind(".dump") != string::npos) {
        read_dump_file();
        return;
    }
    if (trajectory_file_name_.rfind(".atom") != string::npos) {
        read_dump_file();
        return;
    }
    if (trajectory_file_name_.rfind(".xyz") != string::npos) {
        read_xyz_file();
        return;
    }
    if (trajectory_file_name_.rfind("XDATCAR") != string::npos) {
        read_vasp_file();
        return;
    }
    if (trajectory_file_type_ == "abc") {
        read_abc_file();
        return;
    }
    
    size_t period_location = trajectory_file_name_.rfind(".");
    if (period_location != string::npos) {
        cerr << "ERROR: The trajectory file type: ";
        cerr << "\033[1;34m";
        cerr << trajectory_file_name_.substr(period_location);
        cerr << "\033[0m";
        cerr << " does not match any known type\n";
        cerr << "     : Please provide .trr, .xtc, .dump, XDATCAR";
        cerr << endl;
        exit(1);
    }
}


// Function to wrap the system trajectory
// if the trajectory provided is unwrapped
void Trajectory::wrap_coordinates()
{
#pragma omp parallel for
    for (std::size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        for (std::size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (std::size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_[i_frame][i_atom][i_dimension] -= box_length_[i_frame][i_dimension] * round( trajectory_[i_frame][i_atom][i_dimension]/ box_length_[i_frame][i_dimension] );
            }
        }
    }
}


// Function to unwrap the system trajectory
// if the trajectory provided is not already unwrapped
void Trajectory::unwrap_coordinates()
{
    vector< vector<double> > tmp_frame (number_of_system_atoms_, vector<double>(dimension_, 0.0));
    double delta_r = 0.0;
    
#pragma omp parallel for
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            tmp_frame[i_atom][i_dimension] = trajectory_[0][i_atom][i_dimension];
        }
    }
    
#pragma omp parallel
{
    for (size_t i_frame = 1; i_frame < end_frame_ - start_frame_; ++i_frame) {
#pragma omp for private(delta_r)
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                delta_r = trajectory_[i_frame][i_atom][i_dimension] - tmp_frame[i_atom][i_dimension];
                delta_r = delta_r - box_length_[i_frame][i_dimension]*round(delta_r/box_length_[i_frame][i_dimension]);
                tmp_frame[i_atom][i_dimension] = trajectory_[i_frame][i_atom][i_dimension];
                trajectory_[i_frame][i_atom][i_dimension] = trajectory_[i_frame-1][i_atom][i_dimension] + delta_r;
            }
        }
#pragma omp barrier
    }
}
}


// Function to compute velocities from a coordinate trajectory using
// finite difference of first order.  This is needed to do VACF for VASP
void Trajectory::compute_velocities()
{
    if (trajectory_delta_time_ == 0.0) {
        trajectory_delta_time_ = 1.0;
    }
    
#pragma omp parallel
{
    double delta_r = 0.0;
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_ - 1; ++i_frame) {
#pragma omp for
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                delta_r = trajectory_[i_frame+1][i_atom][i_dimension] - trajectory_[i_frame][i_atom][i_dimension];
                delta_r = delta_r - box_length_[i_frame][i_dimension]*round(delta_r/box_length_[i_frame][i_dimension]);
                trajectory_[i_frame][i_atom][i_dimension] = delta_r/trajectory_delta_time_;
            }
        }
#pragma omp barrier
    }
}
}



// Function to determine the indexes of
// of the atom type desired
void Trajectory::select_atoms(vector< unsigned int > & selected_atom_indexes, string const & atom_type_to_select, string const & atom_group_to_select)
{
    //check if selected_atom_indexes are empty, otherwise using push_back later may cause problem
    if (!selected_atom_indexes.empty()) {
        selected_atom_indexes.clear();
    }
    for (unsigned int i_atom = 0; i_atom < system_atom_types_.size(); ++i_atom) {
        // compare atom types using strings of equal length
        // this helps distringuish atoms at different levels, e.g., H vs HW vs HW1,HW2
        if (system_atom_types_[i_atom].find(atom_type_to_select) != string::npos || atom_type_to_select == "all") {
            if (atom_group_to_select == "system") {
                selected_atom_indexes.push_back(i_atom);
                continue;
            }
            else if (atom_group_to_select == "solvent" || atom_group_to_select == "SOL") {
                if (molecule_id_[i_atom].find("SOL") != string::npos) {    //TODO: Note that this selection only works for gromacs files
                    selected_atom_indexes.push_back(i_atom);
                }
                continue;
            }
            else if (atom_group_to_select == "non-solvent") {
                if (molecule_id_[i_atom].find("SOL") == string::npos) {
                    selected_atom_indexes.push_back(i_atom);
                }
                continue;
            }
            // add more atom_group types OR add more selection features here
            else {
                cerr << "\n";
                cerr << "ERROR: Unrecognized atom group \"" << atom_group_to_select << "\"" << endl;
                cerr << "     : Options for atom group: system, solvent, non-solvent" << endl;
                exit(1);
            }
        }
    }
    if (selected_atom_indexes.empty()) {
        cerr << "\n";
        cerr << "ERROR: Unrecognized atom type \"" << atom_type_to_select << "\"" << endl;
        cerr << "     : Atom types are provided in GROMACS .gro file (3rd column, e.g., OW, HW) \n";
        cerr << "     : or in LAMMPS trajectory file (e.g., 1, 2) \n" << endl;
        exit(1);
    }
    
    cout << "\nSelected atoms: " << endl;
    cout << atom_type_to_select << " in " << atom_group_to_select << ", " << selected_atom_indexes.size() << endl;
    cout << endl;
}


// Function to read .gro file,
// function stores the molecule ids, atom types
// the trajectory must be read prior to gro file read
#ifdef GROMACS
void Trajectory::read_gro_file()
{
    ifstream gro_file(gro_file_name_, ifstream::in);
    if (!gro_file) {
        cerr << "ERROR: Failed to open .gro file\n" << endl;
        exit(1);
    }
    
    // allocate memory for molecule_id_, system_atom_types_,
    molecule_id_.resize(number_of_system_atoms_);
    system_atom_types_.resize(number_of_system_atoms_);
    
    
    string line;
    getline(gro_file, line); 	// discard first two lines: 'system name' and 'number of atoms'
    getline(gro_file, line);
        
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_ && !gro_file.eof(); ++i_atom) {     
        char * read_id   = new char [10];
        char * read_type = new char [5];
        
        gro_file.read(read_id, 10);
        molecule_id_[i_atom] = read_id;

        gro_file.read(read_type, 5);
        system_atom_types_[i_atom] = read_type;

        getline(gro_file, line);
    }
    
    gro_file.close();
}
#endif


// Function to read trajectory from a .xtc file
#ifdef GROMACS
void Trajectory::read_xtc_file()
{
	// check dimension, xdrfile library only works for 3d
	// DIM = 3 defined in "xdrfile.h"
	if (dimension_ != 3) {
		cerr << "ERROR: .trr file format only works for 3d data\n" << endl;
		exit(1);
	}
    // .xtc file only contains coordinates
    if (trajectory_data_type_ != "coordinate") {
        cerr << "ERROR: .xtc file only contains coordinates, no velocities or forces\n" << endl;
        exit(1);
    }

	// assign number of atoms
	int tmp_number_of_system_atoms = 0;
	read_xtc_natoms(strdup(trajectory_file_name_.c_str()), &tmp_number_of_system_atoms);
	number_of_system_atoms_ = static_cast< unsigned int >(tmp_number_of_system_atoms);
    
    // resize _trajectory, box_length_, average_box_length_
    allocate_memory();
    
    int   step_number = 0;
    float time        = 0.0;
    float output_precision   = 0.0;
    matrix box        = {{0.0}};
    rvec * coordinate = new rvec[tmp_number_of_system_atoms];
    
    XDRFILE * trajectory_file = xdrfile_open(trajectory_file_name_.c_str(), "r");
    if (trajectory_file == NULL) {
        cerr << "ERROR: Failed to open .xtc file\n" << endl;
        exit(1);
    }
    
    // discard frames before Trajectory::start_frame_
    for (size_t i_frame = 0; i_frame < start_frame_; ++i_frame) {
        bool is_file_read = !read_xtc(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, box, coordinate, &output_precision);
        
        if (!is_file_read) {
            cerr << "\n";
            cerr << "ERROR: Failed to read trajectory" << endl;
            cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
            exit(1);
        }
        if (i_frame == 1) {
            trajectory_delta_time_ = static_cast< double >(time);
        }
    }
    
    if (end_frame_ == start_frame_) {
        cerr << "WARNING: trajectory delta time not known since only one frame is read\n" << endl;
    }
    
    // report reading progress every 5 percent read
    unsigned int frame_gap_to_report = static_cast<unsigned int>(0.05 * (end_frame_ - start_frame_));
    if (frame_gap_to_report < 1) {
        frame_gap_to_report = 1;
    }
    
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        bool is_file_read = !read_xtc(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, box, coordinate, &output_precision);
        
        if (!is_file_read) {
            cerr << "\n";
            cerr << "ERROR: Failed to read trajectory" << endl;
            cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
            exit(1);
        }
        
        save_frame(i_frame, time, box, coordinate);
        
        if (is_run_mode_verbose_) {
            if ((i_frame + 1) % frame_gap_to_report == 0) {
                cout << "\r" << (i_frame + 1) << " frames read";
                cout.flush();
            }
        }
    }
    cout << "\r" << (end_frame_ - start_frame_) << " frames read";
    cout.flush();
    cout << endl;
    
    // normalize averaged box lengths
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        average_box_length_[i_dimension] /= end_frame_ - start_frame_;
    }
    
    xdrfile_close(trajectory_file);

    delete [] coordinate;
}
#endif


// Function to read trajectory from a .trr file
#ifdef GROMACS
void Trajectory::read_trr_file()
{
	// check dimension, xdrfile library only works for 3d,
	// DIM = 3 defined in "xdrfile.h"
	if (dimension_ != 3) {
		cerr << "ERROR: .trr file format only works for 3d data\n" << endl;
		exit(1);
	}
	
	// assign number of atoms
	int tmp_number_of_system_atoms;
	read_trr_natoms(strdup(trajectory_file_name_.c_str()), &tmp_number_of_system_atoms);
	number_of_system_atoms_ = static_cast< unsigned int >(tmp_number_of_system_atoms);
    
    // resize _trajectory, box_length_, average_box_length_
    allocate_memory();
    
    int   step_number = 0;
    float time        = 0.0;
    float lambda      = 0.0;
    matrix box        = {{0.0}};
    
    XDRFILE * trajectory_file = xdrfile_open(trajectory_file_name_.c_str(), "r");
    if (trajectory_file == NULL) {
        cerr << "ERROR: Failed to open .trr file\n" << endl;
        exit(1);
    }
    
	// discard frames before Trajectory::start_frame_
	for (size_t i_frame = 0; i_frame < start_frame_; ++i_frame) {
		bool is_file_read = !read_trr(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, &lambda, box, NULL, NULL, NULL);
        
        if (!is_file_read) {
            cerr << "\n";
			cerr << "ERROR: Failed to read trajectory" << endl;
			cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
			exit(1);
		}
        if (i_frame == 1) {
            trajectory_delta_time_ = static_cast< double >(time);
        }
	}

	if (end_frame_ == start_frame_) {
		cerr << "WARNING: trajectory delta time not known since only one frame is read\n" << endl;
	}
	
    // report reading progress every 5 percent read
    unsigned int frame_gap_to_report = static_cast<unsigned int>(0.05 * (end_frame_ - start_frame_));
    if (frame_gap_to_report < 1) {
        frame_gap_to_report = 1;
    }
    
	if (trajectory_data_type_ == "coordinate") {
		rvec * coordinate = new rvec[tmp_number_of_system_atoms];

		for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
			bool is_file_read = !read_trr(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, &lambda, box, coordinate, NULL, NULL);
            
            if (!is_file_read) {
                cerr << "\n";
                cerr << "ERROR: Failed to read trajectory" << endl;
                cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
                exit(1);
            }
			
            save_frame(i_frame, time, box, coordinate);
            
            if (is_run_mode_verbose_) {
                if ((i_frame + 1) % frame_gap_to_report == 0) {
                    cout << "\r" << (i_frame + 1) << " frames read";
                    cout.flush();
                }
            }
        }
        cout << "\r" << (end_frame_ - start_frame_) << " frames read";
        cout.flush();
        delete [] coordinate;
	}
	else if (trajectory_data_type_ == "velocity") {
		rvec * velocity = new rvec[tmp_number_of_system_atoms];

        for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
			bool is_file_read = !read_trr(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, &lambda, box, NULL, velocity, NULL);
            
            if (!is_file_read) {
                cerr << "\n";
                cerr << "ERROR: Failed to read trajectory" << endl;
                cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
                exit(1);
            }
            
            save_frame(i_frame, time, box, velocity);
            
            if (is_run_mode_verbose_) {
                if ((i_frame + 1) % frame_gap_to_report == 0) {
                    cout << "\r" << (i_frame + 1) << " frames read";
                    cout.flush();
                }
            }
		}
        cout << "\r" << (end_frame_ - start_frame_) << " frames read";
        cout.flush();
        delete [] velocity;
	}
	else if (trajectory_data_type_ == "force") {
		rvec * force = new rvec[tmp_number_of_system_atoms];

        for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
			bool is_file_read = !read_trr(trajectory_file, tmp_number_of_system_atoms, &step_number, &time, &lambda, box, NULL, NULL, force);

            if (!is_file_read) {
                cerr << "\n";
                cerr << "ERROR: Failed to read trajectory" << endl;
                cerr << "     : Please check the total number of frames in your trajectory first\n" << endl;
                exit(1);
            }
            
            save_frame(i_frame, time, box, force);
            
            if (is_run_mode_verbose_) {
                if ((i_frame + 1) % frame_gap_to_report == 0) {
                    cout << "\r" << (i_frame + 1) << " frames read";
                    cout.flush();
                }
            }
		}
        cout << "\r" << (end_frame_ - start_frame_) << " frames read";
        cout.flush();
        delete [] force;
    }
    else {
        cerr << "\n";
        cerr << "ERROR: Use of unrecognized trajectory data type: ";
        cerr << "\033[1;25m";
        cerr << trajectory_data_type_;
        cerr << "\033[0m";
        cerr << "     : please check input file\n" << endl;
        exit(1);
    }
    cout << endl;
	
    // normalize averaged box lengths
	for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
		average_box_length_[i_dimension] /= end_frame_ - start_frame_;
	}
	
	xdrfile_close(trajectory_file);
}
#endif


// Function to read trajectory from a XDATCAR file
void Trajectory::read_vasp_file()
{
    if (dimension_ != 3) {
        cerr << "ERROR: VASP trajectory only has three dimensions\n" << endl;
        exit(1);
    }
    
    if (trajectory_data_type_ != "coordinate") {
        cerr << "ERROR: VASP trajectory only contains coordinates.";
        exit(1);
    }
    
    ifstream trajectory_file(trajectory_file_name_);
    if (!trajectory_file) {
        cerr << "ERROR: Failed to open trajectory file: ";
        cerr << "\033[1;25m";
        cerr << trajectory_file_name_;
        cerr << "\033[0m\n";
        cerr << endl;
        exit(1);
    }
    
    string read_word;
    
    getline(trajectory_file, read_word);
    getline(trajectory_file, read_word);
    
    double scaling_factor = stod(read_word);
    
    vector< double > temp_average_box_length (dimension_, 0.0);
    
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        for (size_t lines_to_read = 0; lines_to_read <= i_dimension; ++lines_to_read) {
            trajectory_file >> read_word;
        }
        temp_average_box_length[i_dimension] = stod(read_word)*scaling_factor;
        getline(trajectory_file, read_word);
    }
    
    getline(trajectory_file, read_word);

    {
        stringstream line_stream(read_word);
        vector<string> atom_types_in_file;
        while (line_stream >> read_word) {
            atom_types_in_file.push_back(read_word);
        }
    
        getline(trajectory_file, read_word);

        vector< int > numberofTypesInFile;
        stringstream line_stream2(read_word);
        while (line_stream2 >> read_word) {
            number_of_atoms_of_each_type_.push_back(stoi(read_word));
            number_of_system_atoms_ += stoi(read_word);
        }
    
        system_atom_types_.resize(0);
        
        for (size_t atom_type = 0; atom_type < atom_types_in_file.size(); ++atom_type) {
            for (size_t i_atom = 0; i_atom < number_of_atoms_of_each_type_[atom_type]; ++i_atom) {
                system_atom_types_.push_back(atom_types_in_file[atom_type]);
            }
        }
    }
    
    // resize trajectory_, box_length_, average_box_length_
    allocate_memory();
    
    average_box_length_ = temp_average_box_length;
    
#pragma omp parallel for
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            box_length_[i_frame][i_dimension] = average_box_length_[i_dimension];
        }
    }
    
    // ignore frames from 0 to start_frame_
    for (size_t i_frame = 0; i_frame < start_frame_; ++i_frame) {
        getline(trajectory_file, read_word);
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            getline(trajectory_file, read_word);
        }
    }
    
    // report reading progress every 5 percent read
    unsigned int frame_gap_to_report = static_cast<unsigned int>(0.05 * (end_frame_ - start_frame_));
    if (frame_gap_to_report < 1) {
        frame_gap_to_report = 1;
    }

    bool is_direct_coord = false;
    
    // read frames from start_frame_ to end_frame_
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        getline(trajectory_file, read_word);
        is_direct_coord = false;
        if (read_word[0] == 'D' || read_word[0] == 'd') {
            is_direct_coord = true;
        }
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file >> read_word;
                if (is_direct_coord) {
                    trajectory_[i_frame][i_atom][i_dimension] = stod(read_word)*average_box_length_[i_dimension];
                }
                else {
                    trajectory_[i_frame][i_atom][i_dimension] = stod(read_word);
                }
            }
        }
        getline(trajectory_file, read_word);
        
        if (is_run_mode_verbose_) {
            if ((i_frame + 1) % frame_gap_to_report == 0) {
                cout << "\r" << (i_frame + 1) << " frames read";
                cout.flush();
            }
        }
    }
    cout << "\r" << (end_frame_ - start_frame_) << " frames read";
    cout.flush();
    cout << endl;
    
    trajectory_file.close();
}


// Function to read trajectory from a .dump file
void Trajectory::read_dump_file()
{
    if (dimension_ != 3) {
        cerr << "ERROR: LAMMPS trajectory dump must have 3 dimensions.\n" << endl;
        exit(1);
    }
    ifstream trajectory_file(trajectory_file_name_);
    
    vector< vector< double > > box_dimensions(dimension_, vector< double > (2, 0.0));
    
    if (!trajectory_file) {
        cerr << "ERROR: LAMMPS dump trajectory File not found\n" << endl;
        exit(1);
    }

    unsigned int current_frame_time_step;
    bool is_scaled = false;

    for (size_t i_frame = 0; i_frame < end_frame_; ++i_frame) {
        // ----------------------- Start of dump file read ------------------------------
        //      :: Add checks for wrapped trajectory etc. (atom type, xyz etc.)
        
        string line_stream;
        getline(trajectory_file, line_stream); // Skip first line with ITEM: TIMESTEP
        trajectory_file >> current_frame_time_step;
        getline(trajectory_file, line_stream);
        
        getline(trajectory_file, line_stream); // Skip line with ITEM: NUMBER OF ATOMS
        trajectory_file >> number_of_system_atoms_;
        getline(trajectory_file, line_stream);
        
        getline(trajectory_file, line_stream); // SKIP line with ITEM: BOX BOUNDS ...
        
        if (i_frame == start_frame_) {
            // resize _trajectory, box_length_, average_box_length_
            allocate_memory();
        }
        
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_file >> box_dimensions[i_dimension][0] >> box_dimensions[i_dimension][1];
            double dimension_length = box_dimensions[i_dimension][1] - box_dimensions[i_dimension][0];
            average_box_length_[i_dimension] += dimension_length;
            box_length_[i_frame][i_dimension] = dimension_length;
        }
        getline(trajectory_file, line_stream);
        getline(trajectory_file, line_stream);
        

        if (i_frame == 0) {
            if (trajectory_data_type_ == "coordinate") {
                if (line_stream.find("xu yu zu") == string::npos) { // Check if dump has right format
                    if (line_stream.find("x y z") == string::npos) {
                        if (line_stream.find("xs ys zs") == string::npos) {
                            if (line_stream.find("xsu ysu zsu") == string::npos) {
                                cerr << "Unknown coordinates(x,y,z) format. \n";
                                cerr << "Trajectory file must contain one of the following formats. \n";
                                cerr << "id type xs ys zs" << endl;
                                cerr << "id type xsu ysu zsu" << endl;
                                cerr << "id type x y z" << endl;
                                cerr << "id type xu yu zu" << endl;
                                exit(1);
                            }
                            is_scaled = true;
                        }
                        is_scaled = true;
                    }
                }
            }
            else if (trajectory_data_type_ == "velocity") {
                size_t format_check = line_stream.find("vx vy vz");
                if (format_check == string::npos) { // Check if dump has right format, with velocity
                    cerr << "Incorrect velocity(vx,vy,vz) format. \n";
                    cerr << "Dump file does not contain atom velocities. \n";
                    cerr << "Correct format :: id type vx vy vz" << endl;
                    exit(1);
                }
            }
            else if (trajectory_data_type_ == "force") {
                size_t format_check = line_stream.find("fx fy fz");
                if (format_check == string::npos) { // Check if dump has right format, with force
                    cerr << "Incorrect force(fx,fy,fz) format. \n";
                    cerr << "Dump file does not contain forces on atoms. \n";
                    cerr << "Correct format :: id type fx fy fz" << endl;
                    exit(1);
                }
            }
            else {
                cerr << "ERROR: Unrecognized LAMMPS trajectory format: " << line_stream << "\n";
                cerr << "     : Dump file must specify coordinate/velocity/force type." << endl;
                exit(1);
            }
        }
        
        // ignore frames from 0 to start_frame_
        if (i_frame < start_frame_) {
            for (size_t atom_index = 0; atom_index < number_of_system_atoms_; ++atom_index) {
                getline(trajectory_file, line_stream);
            }
            continue;
        }

        system_atom_types_.resize(number_of_system_atoms_); // reinitialize system_atom_types_ vector

        unsigned int atom_number_index = 0;
        unsigned int tmp_atom_type_variable = 0;
        
        if (!is_scaled) {
            for (size_t atom_index = 0; atom_index < number_of_system_atoms_; ++atom_index) {
                if (i_frame == start_frame_) {
                    trajectory_file >> atom_number_index;
                    trajectory_file >> system_atom_types_[atom_number_index - 1];
                }
                else {
                    trajectory_file >> atom_number_index >> tmp_atom_type_variable;
                }
                for (size_t dimension_index = 0; dimension_index < dimension_; ++dimension_index) {
                    trajectory_file >> trajectory_[i_frame - start_frame_][atom_number_index - 1][dimension_index];
                }
            }
        }
        else {
            for (size_t atom_index = 0; atom_index < number_of_system_atoms_; ++atom_index) {
                if (i_frame == start_frame_) {
                    trajectory_file >> atom_number_index;
                    trajectory_file >> system_atom_types_[atom_number_index - 1];
                }
                else {
                    trajectory_file >> atom_number_index >> tmp_atom_type_variable;
                }
                for (size_t dimension_index = 0; dimension_index < dimension_; ++dimension_index) {
                    trajectory_file >> trajectory_[i_frame][atom_number_index - 1][dimension_index];
                    trajectory_[i_frame - start_frame_][atom_number_index - 1][dimension_index] *= box_length_[i_frame][dimension_index];
                }
            }
        }
        
        getline(trajectory_file, line_stream);
        
        if (i_frame % 50 == 0) {
            cout << "\r" << (i_frame + 1) << " frames read";
            cout.flush();
        }
    }
    cout << endl;
    cout.flush();
    
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        average_box_length_[i_dimension] /= end_frame_ - start_frame_;
    }
    
    trajectory_file.close();
}


// Function to read trajectory from a .xyz file
void Trajectory::read_xyz_file()
{
    ifstream trajectory_file(trajectory_file_name_);
    if (!trajectory_file) {
        cerr << "ERROR: Failed to open trajectory file: ";
        cerr << "\033[1;25m";
        cerr << trajectory_file_name_;
        cerr << "\033[0m\n";
        cerr << endl;
        exit(1);
    }
    
    string read_word;
    
    bool is_boxlength_set = false;
    if (average_box_length_.size() > 0) {
        is_boxlength_set = true;
    }
    
    // set number of system atoms
    getline(trajectory_file, read_word);
    number_of_system_atoms_ = stoi(read_word);
    getline(trajectory_file, read_word);
    
    // resize trajectory_, box_length_, average_box_length_
    allocate_memory();
    
    // determine the system atom types and find boxlength if not given
    system_atom_types_.resize(number_of_system_atoms_);
    vector<double> box_low(dimension_, 0.0);
    vector<double> box_high(dimension_, 0.0);
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        trajectory_file >> read_word;
        system_atom_types_[i_atom] = read_word;
        
        // check if this is near the edge of the box
        if (!is_boxlength_set) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file >> read_word;
                if (stod(read_word) < box_low[i_dimension]) {
                    box_low[i_dimension] = stod(read_word);
                }
                if (stod(read_word) > box_high[i_dimension]) {
                    box_high[i_dimension] = stod(read_word);
                }
            }
        }
        
        getline(trajectory_file, read_word);
    }
    trajectory_file.clear();
    trajectory_file.seekg(0);
    
    // set the box length if not already set
    if (!is_boxlength_set) {
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            average_box_length_[i_dimension] = box_high[i_dimension] - box_low[i_dimension];
            cout << average_box_length_[i_dimension] << " ";
        }
        cout << endl;
    }
    
    // report reading progress every 5 percent read
    unsigned int frame_gap_to_report = static_cast<unsigned int>(0.05 * (end_frame_ - start_frame_));
    if (frame_gap_to_report < 1) {
        frame_gap_to_report = 1;
    }
    
    for (size_t i_frame = 0; i_frame < start_frame_; ++i_frame) {
        getline(trajectory_file, read_word);
        getline(trajectory_file, read_word);
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            getline(trajectory_file, read_word);
        }
    }
    
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        // ----------------------- Start of dump file read ------------------------------
        //      :: Add checks for wrapped trajectory etc. (atom type, xyz etc.)
        if (trajectory_file.eof()) {
            cerr << endl;
            cerr << "ERROR: More frames were requested then exists in the trajectory file\n";
            cerr << "     : Please check your input file and trajectory file for consistency\n";
            cerr << endl;
            exit(1);
        }
        getline(trajectory_file, read_word);
        getline(trajectory_file, read_word);
        
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            trajectory_file >> read_word;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file >> read_word;
                trajectory_[i_frame][i_atom][i_dimension] = stod(read_word);
            }
        }
        
        trajectory_file >> read_word;
        
        if (is_boxlength_set) {
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                box_length_[i_frame][i_dimension] = average_box_length_[i_dimension];
            }
        }
        
        if (is_run_mode_verbose_) {
            if ((i_frame + 1) % frame_gap_to_report == 0) {
                cout << "\r" << (i_frame + 1) << " frames read";
                cout.flush();
            }
        }
    }
    cout << "hello" << endl;

    trajectory_file.close();
}


// Function read the location files generated by Nate's abc sampling gromacs mod
void Trajectory::read_abc_file()
{
    if (dimension_ != 3) {
        cerr << "ERROR: ABC trajectory only has three dimensions\n" << endl;
        exit(1);
    }
    
    if (trajectory_data_type_ != "coordinate") {
        cerr << "ERROR: ABC trajectory only contains coordinates.";
        exit(1);
    }
    
    ifstream trajectory_file(trajectory_file_name_);
    if (!trajectory_file) {
        cerr << "ERROR: Failed to open trajectory file: ";
        cerr << "\033[1;25m";
        cerr << trajectory_file_name_;
        cerr << "\033[0m\n";
        cerr << endl;
        exit(1);
    }
    
    string read_word;
    
    // set number of system atoms
    trajectory_file >> read_word;
    number_of_system_atoms_ = stoi(read_word);
    getline(trajectory_file, read_word);
    getline(trajectory_file, read_word);
    getline(trajectory_file, read_word);
    
    // determine the system atom types
    system_atom_types_.resize(number_of_system_atoms_);
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        trajectory_file >> read_word;
        system_atom_types_[i_atom] = read_word;
        getline(trajectory_file, read_word);
    }
    trajectory_file.clear();
    trajectory_file.seekg(0);
    
    // resize trajectory_, box_length_, average_box_length_
    allocate_memory();
    
    // report reading progress every 5 percent read
    unsigned int frame_gap_to_report = static_cast<unsigned int>(0.05 * (end_frame_ - start_frame_));
    if (frame_gap_to_report < 1) {
        frame_gap_to_report = 1;
    }
    
    getline(trajectory_file, read_word);
    
    // Read in the trajectory
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        if (trajectory_file.eof()) {
            cerr << endl;
            cerr << "ERROR: More frames were requested then exists in the trajectory file\n";
            cerr << "     : Please check your input file and trajectory file for consistency\n";
            cerr << endl;
            exit(1);
        }

        getline(trajectory_file, read_word);
        getline(trajectory_file, read_word);
        
        for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
            trajectory_file >> read_word;
            for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
                trajectory_file >> read_word;
                trajectory_[i_frame][i_atom][i_dimension] = stod(read_word);
            }
        }
        
        trajectory_file >> read_word;

        if (is_run_mode_verbose_) {
            if ((i_frame + 1) % frame_gap_to_report == 0) {
                cout << "\r" << (i_frame + 1) << " frames read";
                cout.flush();
            }
        }
        
        // set boxlength equal to average_box_length
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            box_length_[i_frame][i_dimension] = average_box_length_[i_dimension];
        }
    }
    cout << "\r" << (end_frame_ - start_frame_) << " frames read";
    cout << endl;
    
    trajectory_file.close();
}


// Inline function to store the extracted coordinates from a .trr or .xtc file
// into the trajectory_ variable.  This function is needed since the library
// for reading .trr/.xtc requires a typedef rvec and we use a vector for the trajectory data
#ifdef GROMACS
void Trajectory::save_frame(unsigned int const & i_frame, float const & time, matrix const & box, rvec const * coordinate)
{
    // store the time step
    if (i_frame == 1) {
        trajectory_delta_time_ = static_cast<double>(time);
    }
    
    // store the box lengths
    for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
        // assume box is a right rectangular prism
        box_length_[i_frame][i_dimension] = box[i_dimension][i_dimension];
        average_box_length_[i_dimension] += box[i_dimension][i_dimension];
    }
    
    // convert coordinate data to Trajectory::trajectory_
    // Todo: may improve this conversion later
#pragma omp parallel for
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        for (size_t i_dimension = 0; i_dimension < dimension_; ++i_dimension) {
            trajectory_[i_frame][i_atom][i_dimension] = coordinate[i_atom][i_dimension];
        }
    }
}
#endif

// inline function to set the atom types and groups to all and system
// This function is for monatomic systems so that a gro file is not always
// required by the user
#ifdef GROMACS
void Trajectory::set_monatomic()
{
    // allocate memory for molecule_id_, system_atom_types_,
    molecule_id_.resize(number_of_system_atoms_);
    system_atom_types_.resize(number_of_system_atoms_);
    
#pragma omp parallel for
    for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
        molecule_id_[i_atom] = "system";
        system_atom_types_[i_atom] = "all";
    }
}
#endif


// allocate memory for trajectory_, box_length_, and average_box_length_
void Trajectory::allocate_memory()
{
    // allocate memory for trajectory_
    try {
        trajectory_ = new double ** [end_frame_ - start_frame_];
#pragma omp parallel for
        for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
            trajectory_[i_frame] = new double * [number_of_system_atoms_];
            for (size_t i_atom = 0; i_atom < number_of_system_atoms_; ++i_atom) {
                trajectory_[i_frame][i_atom] = new double [dimension_];
            }
        }
    }
    catch (...) {
        cerr << "ERROR: The amount of memory needed to store all of the frames in the trajectory\n";
        cerr << "     : exceeded the amount of memory available to this program.";
        cerr << endl;
        exit(1);
    }
    
    // allocate memory for box_length_
    box_length_.resize(end_frame_ - start_frame_);
    for (size_t i_frame = 0; i_frame < end_frame_ - start_frame_; ++i_frame) {
        box_length_[i_frame].resize(dimension_);
    }
    
    // allocate memory for average_box_length_
    average_box_length_.resize(dimension_, 0.0);
}
