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
#include "StructureFactor_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "StructureFactor.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
	
    StructureFactor structure_factor;
	
    structure_factor.read_command_inputs(argc, argv);
    structure_factor.read_input_file();
    structure_factor.read_trajectory();
    structure_factor.compute_S_k();
    structure_factor.write_S_k();
    
    cout << "Successfully computed structure factor.\n";
    cout << "Cheers. \n";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--              Structure Factor              --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: S_k.in)\n";
}