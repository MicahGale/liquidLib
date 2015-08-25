//
//  MeanSquaredDisplacement.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "MeanSquaredDisplacement_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "MeanSquaredDisplacement.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
    MeanSquaredDisplacement mean_squared_displacement;
    
    mean_squared_displacement.read_command_inputs(argc, argv);
    mean_squared_displacement.read_input_file();
    mean_squared_displacement.read_trajectory();
    mean_squared_displacement.compute_r2_t();
    mean_squared_displacement.write_r2_t();
    
    cout << "Successfully computed mean squared displacement\n";
    cout << "Aloha";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--          Mean Squared Displacement         --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: r2_t.in)\n";
    cout << "\n";
}