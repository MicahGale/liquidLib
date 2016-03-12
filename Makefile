SHELL=/bin/bash

##########################
## USER OPTION SETTINGS ##
##########################
# Set to “yes” if OpenMP or xdrfile library should be used
USE_OMP    = no
USE_XDRLIB = yes
USE_GSL    = no
USE_BOOST  = no

# Locations of the xdrfile library if USE_XDRLIB is yes
ifeq ($(USE_XDRLIB), yes)
	XDRHEAD = ./xdrfileinclude/
	XDRLIB  = ./xdrfilelib/
endif

######################
## END USER OPTIONS ##
######################


#######################
## DEVELOPER OPTIONS ##
#######################
# These sections do not need to be edited, 
# unless you wish to change compiler options

# Set important directories
SRC     = ./src
INCLUDE = ./include

# determine location of all object files
CPPS   = $(wildcard $(SRC)/*.cpp)
OBJS   = $(CPPS:.cpp=.o)
SOURCE = $(CPPS:.cpp=)

# set location of all derived files
EXECUTABLE  = $(wildcard $(SRC)/*_main.cpp)
EXECUTABLE := $(EXECUTABLE:./src/%_main.cpp=%)

# set system compiler and flags
CXX      = g++
CXXFLAGS = -c -Wall -std=c++11 -I$(INCLUDE) -O3

MAINFLAGS =

#########################
## OPTIONAL COMPONENTS ##
#########################

## Include OpenMP headers and flags
ifeq ($(USE_OMP), yes)
	CXXFLAGS  := $(CXXFLAGS) -DOMP -fopenmp
	MAINFLAGS := $(MAINFLAGS) -DOMP -fopenmp
endif
##

## Determine if gromacs version and libraries should be included
ifeq ($(USE_XDRLIB), yes)
	CXXFLAGS := $(CXXFLAGS) -DGROMACS
	#
	# add directory with libxdrfile.a in it, for an example we placed the header files in 
	# “./xdrfileinclude” and the library file in “./xdrfilelib”
	#
	CXXFLAGS  := $(CXXFLAGS) -I$(XDRHEAD)
	MAINFLAGS := $(MAINFLAGS) -lxdrfile -L$(XDRLIB)
endif
##

## Determine if gsl is to be used
ifeq ($(USE_GSL), yes)
	CXXFLAGS := $(CXXFLAGS) -DGSL
	#
	# add directories where the gsl headers and libraries were installed
	#
	CXXFLAGS  := $(CXXFLAGS) -I/opt/local/include
	MAINFLAGS := $(MAINFLAGS) -I/opt/local/include -lgsl -lgslcblas -L/opt/local/lib
endif
##

## Determine if boost is to be used
ifeq ($(USE_BOOST), yes)
	CXXFLAGS := $(CXXFLAGS) -DBOOST
	#
	# add directories where the boost headers and libraries were installed
	#
	CXXFLAGS  := $(CXXFLAGS) -I/opt/local/include
	MAINFLAGS := $(MAINFLAGS) -I/opt/local/include
endif
##

# set base objects
BASE_OBJECTS = ./src/Trajectory.o ./src/MeanSquaredDisplacement.o

PHONY: all
all: MakeBin depends $(OBJS) $(EXECUTABLE)

# Build executable files
$(EXECUTABLE): %: ./src/%_main.o ./src/%.o $(BASE_OBJECTS)
	$(CXX) $^ -o ./bin/compute$@ $(MAINFLAGS)

# Build all object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up all built files
clean:
	rm -rf $(SRC)/*.o ./bin $(SRC)/depends.mk

# clean only the object files
clean-o:
	rm -rf $(SRC)/*.o

# Make binary directory
MakeBin:
	@ mkdir -p ./bin

# Determine dependencies 
depends:
	@ rm -f $(SRC)/depends.mk
	@ for file in $(SOURCE); do \
		$(CXX) $(CXXFLAGS) -MM $$file.cpp -I$(INCLUDE) -MT $$file.o >> $(SRC)/depends.mk; \
	  done

# Include dependencies
-include depends.mk
