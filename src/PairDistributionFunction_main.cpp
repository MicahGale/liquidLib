//
//  PairDistributionFunction_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "PairDistributionFunction_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "PairDistributionFunction.hpp"

using namespace std;

int main(int argc, char * argv[])
{
	print_executable_header();
	
    PairDistributionFunction pair_distribution_function;
    
    pair_distribution_function.read_command_inputs(argc, argv);
	pair_distribution_function.read_input_file();
	pair_distribution_function.read_trajectory();
	pair_distribution_function.compute_g_r();
    pair_distribution_function.write_g_r();
    
    cout << "Successfully computed pair distribution function\n";
    cout << "bai bai la";
    cout << endl;
    
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--         Pair Distribution Function         --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: g_r.in)\n";
    cout << "\n";
}