//
//  CoherentVanHoveFunction_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "CoherentVanHoveFunction_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "CoherentVanHoveFunction.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
	
    CoherentVanHoveFunction coherent_van_hove_function;
    
    coherent_van_hove_function.read_command_inputs(argc, argv);
    coherent_van_hove_function.read_input_file();
    coherent_van_hove_function.read_trajectory();
    coherent_van_hove_function.compute_G_rt();
    coherent_van_hove_function.write_G_rt();
    
    cout << "Successfully computed coherent van hove correlation\n";
    cout << "Peace Out.";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--         Coherent Van Hove Function         --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: G_rt.in)\n";
    cout << "\n";
}