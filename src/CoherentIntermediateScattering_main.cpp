//
//  CoherentIntermediateScattering_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "CoherentIntermediateScattering_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "CoherentIntermediateScattering.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
	
    CoherentIntermediateScattering Coherent_intermediate_scattering;
    
    Coherent_intermediate_scattering.read_command_inputs(argc, argv);
    Coherent_intermediate_scattering.read_input_file();
    Coherent_intermediate_scattering.read_trajectory();
    Coherent_intermediate_scattering.compute_F_kt();
    Coherent_intermediate_scattering.write_F_kt();
    
    cout << "Successfully computed Coherent intermediate scattering function\n";
    cout << "Cya";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--  Coherent Intermediate Scattering Function --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: F_kt.in)\n";
    cout << "\n";
}