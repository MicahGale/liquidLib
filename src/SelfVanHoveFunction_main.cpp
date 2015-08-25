//
//  SelfVanHoveFunction_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "SelfVanHoveFunction_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "SelfVanHoveFunction.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
	
    SelfVanHoveFunction self_van_hove_function;
    
    self_van_hove_function.read_command_inputs(argc, argv);
    self_van_hove_function.read_input_file();
    self_van_hove_function.read_trajectory();
    self_van_hove_function.compute_Gs_rt();
    self_van_hove_function.write_Gs_rt();
    
    cout << "Successfully computed self van hove correlation\n";
    cout << "Peace Out.";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--           Self Van Hove Function           --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: Gs_rt.in)\n";
    cout << "\n";
}