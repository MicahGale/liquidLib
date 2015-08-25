//
//  SelfIntermediateScattering_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "SelfIntermediateScattering_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "SelfIntermediateScattering.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
	
    SelfIntermediateScattering self_intermediate_scattering;
    
    self_intermediate_scattering.read_command_inputs(argc, argv);
    self_intermediate_scattering.read_input_file();
    self_intermediate_scattering.read_trajectory();
    self_intermediate_scattering.compute_Fs_kt();
    self_intermediate_scattering.write_Fs_kt();
    
    cout << "Successfully computed self intermediate scattering function\n";
    cout << "Cya";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--    Self Intermediate Scattering Function   --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: Fs_kt.in)\n";
    cout << "\n";
}