//
//  FourPointCorrelation_main.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "FourPointCorrelation_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "FourPointCorrelation.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
    FourPointCorrelation four_point_correlation;
    
    four_point_correlation.read_command_inputs(argc, argv);
    four_point_correlation.read_input_file();
    four_point_correlation.read_trajectory();
    four_point_correlation.compute_chi4_t();
    four_point_correlation.write_chi4_t();
    
    cout << "Successfully computed four point correlation\n";
    cout << "Toodles";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--       Four Point Correlation Function      --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: chi4_t.in)\n";
    cout << "\n";
}