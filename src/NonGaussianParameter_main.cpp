//
//  NonGaussianParameter_main.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "NonGaussianParameter_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "NonGaussianParameter.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
    NonGaussianParameter non_gaussian_parameter;
    
    non_gaussian_parameter.read_command_inputs(argc, argv);
    non_gaussian_parameter.read_input_file();
    non_gaussian_parameter.set_output_file_name();
    non_gaussian_parameter.read_trajectory();
    non_gaussian_parameter.compute_alpha2_t();
    non_gaussian_parameter.write_alpha2_t();
    
    cout << "Successfully computed non gaussian parameter and mean squared displacement\n";
    cout << "Bye bye.";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--        Mean Squared Displacement and       --\n";
    cout << "--            Non Gaussian Parameter          --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: alpha2_t.in)\n";
    cout << "\n";
}