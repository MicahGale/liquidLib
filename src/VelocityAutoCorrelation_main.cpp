//
//  VelocityAutoCorrelation.cpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "VelocityAutoCorrelation_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "VelocityAutoCorrelation.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
    VelocityAutoCorrelation velocityautocorrelation;
    
    velocityautocorrelation.read_command_inputs(argc, argv);
    velocityautocorrelation.read_input_file();
    velocityautocorrelation.read_trajectory();
    velocityautocorrelation.compute_vacf_t();
    velocityautocorrelation.write_vacf_t();
    
    cout << "Successfully computed velocity autocorrelaiton function\n";
    cout << "Adios.";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--     Velocity AutoCorrelation Function      --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i :input file name (default input file: vacf_t.in)\n";
    cout << "\n";
}