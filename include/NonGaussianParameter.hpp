//
//  NonGaussianParameter.hpp
//  
//  Copyright (c) 2015 Zhang-Group. All rights reserved.
//  This software is distributed under the MIT license
//  -----------------------------------------------------
//  Contributing authors: Zhikun Cai,
//                        Abhishek Jaiswal,
//                        Nathan Walter
//  -----------------------------------------------------
//
#ifndef LiquidLib_NonGaussianParameter_hpp
#define LiquidLib_NonGaussianParameter_hpp

#include <vector>
#include <string>

#include "MeanSquaredDisplacement.hpp"
#include "Trajectory.hpp"

using namespace std;

class NonGaussianParameter : public MeanSquaredDisplacement {
public:
    NonGaussianParameter();
    virtual ~NonGaussianParameter();
    
    void set_output_file_name();
    void compute_alpha2_t();
    void write_alpha2_t();
    
protected:
    
private:
    inline void print_status(size_t & status);
    
//private variables
    vector< double > alpha2_t_;
};

#endif // defined (LiquidLib_NonGaussianParameter_hpp)
