SHELL=/bin/bash

CC=g++
CCFLAGS=-c -Wall -std=c++11 -I./include/ -O3 

MAINFLAGS=

CMD=	Trajectory.cpp \
	MeanSquaredDisplacement.cpp 	MeanSquaredDisplacement_main.cpp \
	NonGaussianParameter.cpp 	NonGaussianParameter_main.cpp \
	PairDistributionFunction.cpp 	PairDistributionFunction_main.cpp \
	SelfVanHoveFunction.cpp 	SelfVanHoveFunction_main.cpp \
	FourPointCorrelation.cpp 	FourPointCorrelation_main.cpp \
	SelfIntermediateScattering.cpp 	SelfIntermediateScattering_main.cpp \
	VelocityAutoCorrelation.cpp 	VelocityAutoCorrelation_main.cpp \
	StructureFactor.cpp 		StructureFactor_main.cpp

EXECUTABLE=	MeanSquaredDisplacement NonGaussianParameter PairDistributionFunction \
		SelfVanHoveFunction 	FourPointCorrelation SelfIntermediateScattering \
		VelocityAutoCorrelation StructureFactor

INCLUDE=$(addprefix ./include/, $(CMD))
SRC=$(addprefix ./src/, $(CMD))
OBJECTS=$(CMD:.cpp=.o)
BASE_OBJECTS=./src/Trajectory.o ./src/MeanSquaredDisplacement.o

#########################
## OPTIONAL COMPONENTS ##
#########################

## Uncomment following line for parallelized version
#
#CCFLAGS:=$(CCFLAGS) -DOMP -fopenmp
#MAINFLAGS:=$(MAINFLAGS) -DOMP -fopenmp
#
##

## Uncomment following lines for gromacs version
#
CCFLAGS:=$(CCFLAGS) -DGROMACS
#
# add directory with libxdrfile.a in it, for an example the following lines are what we
# used for a directory structure
#
CCFLAGS:=$(CCFLAGS) -I./xdrfileinclude/
MAINFLAGS:=$(MAINFLAGS) -lxdrfile -L./xdrfilelib/
#
##

PHONY: all
all:MakeBin $(OBJECTS) $(EXECUTABLE)

MakeBin:
	mkdir -p ./bin

$(EXECUTABLE): %: ./src/%_main.o ./src/%.o $(BASE_OBJECTS)
	$(CC) $^ -o ./bin/compute$@ $(MAINFLAGS) 

$(OBJECTS): %.o: ./src/%.cpp
	$(CC) $(CCFLAGS) $< -o ./src/$@

clean:
	rm -rf */*.o ./bin

clean-o:
	rm -rf */*.o
