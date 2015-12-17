//
//  BondOrderParameter_main.cpp
//
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#include "BondOrderParameter_main.hpp"

#include <iostream>
#include <string>
#include <cstring>

#include "BondOrderParameter.hpp"

using namespace std;

int main(int argc, char * argv[])
{
    print_executable_header();
    
#if defined(GSL) || defined(BOOST)
    BondOrderParameter bond_order_parameter;
    
    bond_order_parameter.read_command_inputs(argc, argv);
    bond_order_parameter.read_input_file();
    bond_order_parameter.read_trajectory();
    bond_order_parameter.compute_BOP();
    bond_order_parameter.write_BOP();
    
    cout << "Successfully computed Bond Order Parameter\n";
    cout << "Au revoir";
    cout << endl;
#else
    cerr << "ERROR: Bonded Order Parameter can only be computed with\n";
    cerr << "     : Gnu Scientific Library or BOOST Installed.  Please edit\n";
    cerr << "     : Makefile to account for this requirement to continue\n";
    cerr << endl;
    exit(1);
#endif
    return 0;
}


void print_executable_header()
{
    cout << "------------------------------------------------\n";
    cout << "                   LiquidLib                    \n";
    cout << "------------------------------------------------\n";
    cout << "--            Bond Order Parameter            --\n";
    cout << "------------------------------------------------\n";
    cout << "------------------------------------------------\n";
    cout << "-i: input file name (default input file: BOP.in)\n";
    cout << "\n";
}