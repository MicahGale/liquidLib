# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:01:27 2016

@author: Nathan
"""

import Tkinter     
import os    

quantity_options = ("Pair Distribution Function",
                    "Structure Factor",
                    "Bond Order Parameter",                    
                    "Mean Squared Displacement",
                    "Non-Gaussian Parameter",
                    "Four-Point Correlation",
                    "Velocity Auto Correlation",
                    "Self van Hove Correlation",
                    "Coherent van Hove Correlation",
                    "Self Intermediate Scattering Function",
                    "Coherent Intermediate Scattering Function")    

def option_changed(*args):
    if quantity.get() ==  "Pair Distribution Function":
        hide_computation_widgets()
        show_atom_options()
        show_atom2_options()
        show_bin_options()
        
    elif quantity.get() == "Structure Factor":
        hide_computation_widgets()
        show_atom_options()
        show_scattering_options() 
        
    elif quantity.get() == "Bond Order Parameter":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_bop_options()
        
    elif quantity.get() == "Mean Squared Displacement":
        hide_computation_widgets()
        show_atom_options()
        show_time_options() 
        
    elif quantity.get() == "Non-Gaussian Parameter":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()   
        
    elif quantity.get() == "Four-Point Correlation":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_chi4_options()
        
    elif quantity.get() == "Velocity Auto Correlation":
        hide_computation_widgets()
        show_atom_options()
        show_time_options() 
        
    elif quantity.get() == "Self van Hove Correlation":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_bin_options()  
        
    elif quantity.get() == "Coherent van Hove Correlation":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_bin_options()  
        
    elif quantity.get() == "Self Intermediate Scattering Function":
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_scattering_options()  
        
    elif quantity.get() == "Coherent Intermediate Scattering Function":
        hide_computation_widgets()
        hide_computation_widgets()
        show_atom_options()
        show_time_options()
        show_scattering_options()
        
    else:
        return

    return
    
def show_atom_options():
    atomtypes_frame.pack()
    atomtypes_label.pack( side = Tkinter.LEFT)
    atomtypes_value.configure(width = 40)
    atomtypes_value.pack( side = Tkinter.LEFT)
    
    scatteringlengths_frame.pack()
    scatteringlengths_label.pack( side = Tkinter.LEFT)
    scatteringlengths_value.configure(width = 40)
    scatteringlengths_value.pack( side = Tkinter.LEFT)
    
    atomgroup_frame.pack()    
    atomgroup_label.pack( side = Tkinter.LEFT)
    atomgroup_value.configure(width = 40)
    atomgroup_value.pack( side = Tkinter.LEFT)    
    
def show_atom2_options():
    atomtypes2_frame.pack()
    atomtypes2_label.pack( side = Tkinter.LEFT)
    atomtypes2_value.configure(width = 40)
    atomtypes2_value.pack( side = Tkinter.LEFT)
    
    scatteringlengths2_frame.pack()
    scatteringlengths2_label.pack( side = Tkinter.LEFT)
    scatteringlengths2_value.configure(width = 40)
    scatteringlengths2_value.pack( side = Tkinter.LEFT)
    
    atomgroup2_frame.pack()
    atomgroup2_label.pack( side = Tkinter.LEFT)
    atomgroup2_value.configure(width = 40)
    atomgroup2_value.pack( side = Tkinter.LEFT)    

def show_time_options():
    timepoints_frame.pack()
    timepoints_label.pack(side = Tkinter.LEFT)
    timepoints_value.configure(width = 39)
    timepoints_value.pack(side = Tkinter.LEFT)

    timescale_frame.pack()
    timescale_label.pack(side = Tkinter.LEFT)
    timescale_value.configure(width = 39)
    timescale_value.pack(side = Tkinter.LEFT)

    interval_frame.pack()
    interval_label.pack(side = Tkinter.LEFT)
    interval_value.configure(width = 40)
    interval_value.pack(side = Tkinter.LEFT)    
    
def show_bin_options():
    numberbins_frame.pack()
    numberbins_label.pack( side = Tkinter.LEFT)
    numberbins_value.configure(width = 39)
    numberbins_value.pack( side = Tkinter.LEFT)
    
    cutofflength_frame.pack()
    cutofflength_label.pack( side = Tkinter.LEFT)
    cutofflength_value.configure(width = 40)
    cutofflength_value.pack( side = Tkinter.LEFT)    
    
def show_scattering_options():
    kstart_frame.pack()
    kstart_label.pack(side = Tkinter.LEFT)
    kstart_value.configure(width = 39)
    kstart_value.pack(side = Tkinter.LEFT)    
    
    numberbins_frame.pack()
    numberbins_label.pack( side = Tkinter.LEFT)
    numberbins_value.configure(width = 39)
    numberbins_value.pack( side = Tkinter.LEFT)    

## for bop    
def show_bop_options():
    bondorder_frame.pack()
    bondorder_label.pack(side = Tkinter.LEFT)
    bondorder_value.configure(width = 39)
    bondorder_value.pack(side = Tkinter.LEFT)

    cutofflength_frame.pack()
    cutofflength_label.pack( side = Tkinter.LEFT)
    cutofflength_value.configure(width = 40)
    cutofflength_value.pack( side = Tkinter.LEFT)

    qaverage_frame.pack()
    qaverage_label.pack(side = Tkinter.LEFT)
    qaverage_value.configure(width = 40)
    qaverage_value.pack(side = Tkinter.LEFT)

    qbar_frame.pack()
    qbar_label.pack(side = Tkinter.LEFT)
    qbar_value.configure(width = 40)
    qbar_value.pack(side = Tkinter.LEFT)    
    
def show_chi4_options():
    overlap_frame.pack()
    overlap_label.pack(side = Tkinter.LEFT)
    overlap_value.configure(width = 40)
    overlap_value.pack(side = Tkinter.LEFT)

def hide_trajectory_widgets():
    trajectory_frame.pack_forget()
    trajectory_label.pack_forget() 
    trajectory_value.pack_forget()
    
    wrapped_frame.pack_forget()
    wrapped_label.pack_forget()
    wrapped_value.pack_forget()
    
    startframe_frame.pack_forget()
    startframe_label.pack_forget()
    startframe_value.pack_forget()
    
    endframe_frame.pack_forget()
    endframe_label.pack_forget()
    endframe_value.pack_forget()
    
    frameaverage_frame.pack_forget()
    frameaverage_label.pack_forget()
    frameaverage_value.pack_forget()
    
    dimension_frame.pack_forget()
    dimension_label.pack_forget()
    dimension_value.pack_forget()
    
    timestep_frame.pack_forget()
    timestep_label.pack_forget()
    timestep_value.pack_forget()
    
    gro_frame.pack_forget()
    gro_label.pack_forget()
    gro_value.pack_forget()
    
def hide_quantity_widgets():
    quantity_frame.pack_forget()
    quantity_label.pack_forget()
    quantity_dropdown.pack_forget()
    
    verbose_frame.pack_forget()
    verbose_label.pack_forget()
    verbose_value.pack_forget()
    
    outputfile_frame.pack_forget()
    outputfile_label.pack_forget()
    outputfile_value.pack_forget()    
    
def hide_computation_widgets():
    atomtypes_frame.pack_forget()
    atomtypes_label.pack_forget()
    atomtypes_value.pack_forget()

    atomtypes2_frame.pack_forget()
    atomtypes2_label.pack_forget()
    atomtypes2_value.pack_forget()

    scatteringlengths_frame.pack_forget()
    scatteringlengths_label.pack_forget()
    scatteringlengths_value.pack_forget()

    scatteringlengths2_frame.pack_forget()
    scatteringlengths2_label.pack_forget()
    scatteringlengths2_value.pack_forget()

    atomgroup_frame.pack_forget()
    atomgroup_label.pack_forget()
    atomgroup_value.pack_forget()

    atomgroup2_frame.pack_forget()
    atomgroup2_label.pack_forget()
    atomgroup2_value.pack_forget()

    numberbins_frame.pack_forget()
    numberbins_label.pack_forget()
    numberbins_value.pack_forget()

    cutofflength_frame.pack_forget()
    cutofflength_label.pack_forget()
    cutofflength_value.pack_forget()

    kstart_frame.pack_forget()
    kstart_label.pack_forget()
    kstart_value.pack_forget()

    timepoints_frame.pack_forget()
    timepoints_label.pack_forget()
    timepoints_value.pack_forget()

    timescale_frame.pack_forget()
    timescale_label.pack_forget()
    timescale_value.pack_forget()

    interval_frame.pack_forget()
    interval_label.pack_forget()
    interval_value.pack_forget()

    overlap_frame.pack_forget()
    overlap_label.pack_forget()
    overlap_value.pack_forget()

    qaverage_frame.pack_forget()
    qaverage_label.pack_forget()
    qaverage_value.pack_forget()

    qbar_frame.pack_forget()
    qbar_label.pack_forget()
    qbar_value.pack_forget()

    bondorder_frame.pack_forget()
    bondorder_label.pack_forget()
    bondorder_value.pack_forget()
        
def show_quantity_options():
    hide_trajectory_widgets()
    
    quantity_frame.pack()
    quantity_label.pack( side = Tkinter.LEFT )
    quantity_dropdown.configure(width = 39)
    quantity_dropdown.pack()
    
    verbose_frame.pack()
    verbose_label.pack(side = Tkinter.LEFT)
    verbose_value.configure(width = 40)
    verbose_value.pack(side = Tkinter.LEFT)
    
    outputfile_frame.pack(side = Tkinter.BOTTOM, pady = 20)
    outputfile_label.pack(side = Tkinter.LEFT)
    outputfile_value.configure(width = 40)
    outputfile_value.pack(side = Tkinter.LEFT)
    
    option_changed()
        
def show_trajectory_options():
    hide_quantity_widgets()
    hide_computation_widgets()
    
    trajectory_frame.pack()
    trajectory_label.pack( side = Tkinter.LEFT)
    trajectory_value.configure(width = 40)
    trajectory_value.pack( side = Tkinter.LEFT)

    wrapped_frame.pack()
    wrapped_label.pack( side = Tkinter.LEFT)
    wrapped_value.configure(width = 40)
    wrapped_value.pack( side = Tkinter.LEFT)
        
    gro_frame.pack()
    gro_label.pack( side = Tkinter.LEFT)
    gro_value.configure(width = 40)
    gro_value.pack( side = Tkinter.LEFT)   
    
    startframe_frame.pack()
    startframe_label.pack(side = Tkinter.LEFT)
    startframe_value.configure(width = 39)
    startframe_value.pack(side= Tkinter.LEFT) 
    
    endframe_frame.pack()
    endframe_label.pack(side = Tkinter.LEFT)
    endframe_value.configure(width = 39)
    endframe_value.pack(side= Tkinter.LEFT) 
    
    frameaverage_frame.pack()
    frameaverage_label.pack(side = Tkinter.LEFT)
    frameaverage_value.configure(width = 39)
    frameaverage_value.pack(side= Tkinter.LEFT)
    
    timestep_frame.pack()
    timestep_label.pack(side = Tkinter.LEFT)
    timestep_value.configure(width = 40)
    timestep_value.pack(side = Tkinter.LEFT)
    
    dimension_frame.pack()
    dimension_label.pack(side = Tkinter.LEFT)
    dimension_value.configure(width=39)
    dimension_value.pack(side = Tkinter.LEFT)
    
def generate_input_file():
    input_file = file
    if quantity.get() ==  "Pair Distribution Function":
        input_file = open("gofr.in", 'w') 
        write_trajectory_options(input_file)
        write_2atom_options(input_file)
        write_bin_options(input_file)
        
    elif quantity.get() == "Structure Factor":
        input_file = open("sofk.in", 'w')   
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_scattering_options(input_file)
        
    elif quantity.get() == "Bond Order Parameter":
        input_file = open("bop.in", 'w')   
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_bop_options(input_file)
        
    elif quantity.get() == "Mean Squared Displacement":
        input_file = open("msd.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        
    elif quantity.get() == "Non-Gaussian Parameter":
        input_file = open("alpha2.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)        
        
    elif quantity.get() == "Four-Point Correlation":
        input_file = open("chi4.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        write_chi4_options(input_file)
        
    elif quantity.get() == "Velocity Auto Correlation":
        input_file = open("vacf.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        
    elif quantity.get() == "Self van Hove Correlation":
        input_file = open("gsrt.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        write_bin_options(input_file)
        
    elif quantity.get() == "Coherent van Hove Correlation":
        input_file = open("grt.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        write_bin_options(input_file)
        
    elif quantity.get() == "Self Intermediate Scattering Function":
        input_file = open("fskt.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        write_scattering_options(input_file)

    elif quantity.get() == "Coherent Intermediate Scattering Function":
        input_file = open("fkt.in", 'w')        
        write_trajectory_options(input_file)
        write_atom_options(input_file)
        write_time_options(input_file)
        write_scattering_options(input_file)

    else:
        return
            
    input_file.close()
    return
    
def write_trajectory_options(input_file):
    input_file.write("# This file generated by LiquidLib\n")
    input_file.write("# Trajectory Quantities\n")
    input_file.write("is_run_mode_verbose    = " + str(verbose.get()) + "\n")
    input_file.write("trajectory_file_name   = " + str(trajectory_value.get()) + "\n")
    if gro_value.get() != "":
        input_file.write("gro_file_name          = " + str(gro_value.get()) + "\n")
        
    input_file.write("is_wrapped             = " + str(wrapped.get()) + "\n")
    input_file.write("start_frame            = " + str(startframe_value.get()) + "\n")
    input_file.write("end_frame              = " + str(endframe_value.get()) + "\n")
    if timestep_value.get() != "":
        input_file.write("trajectory_delta_time  = " + str(timestep_value.get()) + "\n")
        
    input_file.write("dimension              = " + str(dimension_value.get()) + "\n")
    if outputfile_value.get() == "":
        input_file.write("output_file_name       = output.out\n")
    else:
        input_file.write("output_file_name       = " + str(outputfile_value.get()) + "\n")
        
    input_file.write("number_of_frames_to_average = " + str(frameaverage_value.get()) + "\n")
    
def write_atom_options(input_file):
    input_file.write("\n# Atom Quantities\n")
    if atomtypes_value.get() == "":
        input_file.write("atom_types            = " + "all" + "\n")
    else:
        input_file.write("atom_types            = " + str(atomtypes_value.get()) + "\n")  
          
    if scatteringlengths_value.get() != "":
        input_file.write("scattering_lengths    = " + str(scatteringlengths_value.get()) + "\n")
        
    if atomgroup_value.get() != "":
        input_file.write("atom_group            = " + str(atomgroup_value.get()) + "\n")
        
def write_2atom_options(input_file):
    input_file.write("\n# Atom Quantities\n")
    if atomtypes_value.get() == "":
        input_file.write("atom_types1           = " + "all" + "\n")
    else:
        input_file.write("atom_types1           = " + str(atomtypes_value.get()) + "\n") 
        
    if scatteringlengths_value.get() != "":
        input_file.write("scattering_lengths1   = " + str(scatteringlengths_value.get()) + "\n") 
        
    if atomgroup_value.get() != "":
        input_file.write("atom_group1           = " + str(atomgroup_value.get()) + "\n")
        
    if atomtypes2_value.get() == "":
        input_file.write("atom_types2           = " + "all" + "\n")
    else:
        input_file.write("atom_types2           = " + str(atomtypes2_value.get()) + "\n") 
           
    if scatteringlengths2_value.get() != "":
        input_file.write("scattering_lengths2   = " + str(scatteringlengths2_value.get()) + "\n")
            
    if atomgroup2_value.get() != "":
        input_file.write("atom_group2           = " + str(atomgroup2_value.get()) + "\n")        

def write_time_options(input_file):
    input_file.write("\n# Time Quantities\n")
    if timepoints_value.get() != "":
        input_file.write("number_of_time_points = " + str(timepoints_value.get()) + " \n")
        
    if timescale.get() != "":
        input_file.write("time_scale_type       = " + str(timescale.get()) + "\n")
        
    if interval_value.get() != "":
        input_file.write("frame_interval        = " + str(interval_value.get()) + "\n")   
        
def write_scattering_options(input_file):
    input_file.write("\n# Scattering Quantities\n")
    if numberbins_value.get() != "":
        input_file.write("number_of_bins        = " + str(numberbins_value.get()) + "\n")
        
    if kstart_value.get() != "":
        input_file.write("k_start_index         = " + str(kstart_value.get()) + "\n")
        
def write_bin_options(input_file):
    input_file.write("\n# Binning Quantities\n")
    if numberbins_value.get() != "":
        input_file.write("number_of_bins        = " + str(numberbins_value.get()) + "\n")

    if cutofflength_value.get() != "":
        input_file.write("max_cutoff_length     = " + str(cutofflength_value.get()) + "\n")
        
def write_bop_options(input_file):
    input_file.write("\n# Bond Order Quantities\n")
    input_file.write("bond_parameter_order  = " + str(bondorder_value.get()) + "\n")
    input_file.write("max_cutoff            = " + str(cutofflength_value.get()) + "\n")
    input_file.write("is_averaged           = " + str(is_averaged.get()) + "\n")
    input_file.write("add_bar               = " + str(add_bar.get()) + "\n")
    
def write_chi4_options(input_file):
    input_file.write("overlap_length        = " + str(overlap_value.get()) + "\n")
       
def compute_quantity():
    generate_input_file()
        
    if quantity.get() ==  "Pair Distribution Function":
        os.system('../bin/computePairDistributionFunction -i gofr.in')
        os.remove('gofr.in')
    
    elif quantity.get() == "Structure Factor":
        os.system('../bin/computeStructureFactor - i sofk.in')
        os.remove('sofk.in')
        
    elif quantity.get() == "Bond Order Parameter":
        os.system('../bin/computeBondOrderParameter - i bop.in')
        os.remove('bop.in')
         
    elif quantity.get() == "Mean Squared Displacement":
        os.system('../bin/computeMeanSquaredDisplacement - i msd.in')
        os.remove('msd.in')
        
    elif quantity.get() == "Non-Gaussian Parameter":        
        os.system('../bin/computeNonGaussianParameter - i alpha2.in')
        os.remove('alpha2.in')
        
    elif quantity.get() == "Four-Point Correlation":
        os.system('../bin/computeFourPointCorrelation - i chi4.in')
        os.remove('chi4.in')
        
    elif quantity.get() == "Velocity Auto Correlation":
        os.system('../bin/computeVelocityAutoCorrelation - i vacf.in')
        os.remove('vacf.in')
        
    elif quantity.get() == "Self van Hove Correlation":
        os.system('../bin/computeSelfVanHoveFunction - i gsrt.in')
        os.remove('gsrt.in')
        
    elif quantity.get() == "Coherent van Hove Correlation":
        os.system('../bin/computeCoherentVanHoveFunction - i grt.in')
        os.remove('grt.in')
        
    elif quantity.get() == "Self Intermediate Scattering Function":
        os.system('../bin/computeSelfIntermediateScattering - i fskt.in')
        os.remove('fskt.in')
        
    elif quantity.get() == "Coherent Intermediate Scattering Function":
        os.system('../bin/computeCoherentIntermediateScattering - i fkt.in')
        os.remove('fkt.in')
        
    else:
        return
    
    top.quit()


top = Tkinter.Tk()
top.attributes("-topmost", True)

top.wm_title("LiquidLib")

photo = Tkinter.PhotoImage(file="./LiquidLib_logo.gif")
logo = Tkinter.Label(top, image=photo)
logo.pack()

### Top menu
categories_frame = Tkinter.Frame(top)
categories_frame.pack(pady = 20)
trajectory_option = Tkinter.Button(categories_frame, text = "Trajectory", command = show_trajectory_options)
trajectory_option.pack(in_=categories_frame, side=Tkinter.LEFT, padx = 10)
quantity_option = Tkinter.Button(categories_frame, text = "Quantity", command = show_quantity_options)
quantity_option.pack(in_=categories_frame, side=Tkinter.LEFT, padx = 10)

### Quantity Options
quantity = Tkinter.StringVar(top)
quantity.set("QUANTITY")
quantity.trace("w", option_changed)
quantity_frame = Tkinter.Frame(top)
quantity_label = Tkinter.Label(quantity_frame, text = "Quantity", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
quantity_dropdown = Tkinter.OptionMenu(quantity_frame, quantity, *quantity_options)

verbose = Tkinter.IntVar()
verbose_frame = Tkinter.Frame(top)
verbose_label = Tkinter.Label(verbose_frame, text="", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
verbose_value = Tkinter.Checkbutton(verbose_frame, text = "verbose", justify = Tkinter.LEFT, variable = verbose)

outputfile_frame = Tkinter.Frame(top)
outputfile_label = Tkinter.Label(outputfile_frame, text = "Output File", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
outputfile_value = Tkinter.Entry(outputfile_frame, bd = 5)

### Trajectory Options
trajectory_frame = Tkinter.Frame(top)
trajectory_label = Tkinter.Label(trajectory_frame, text = "Trajectory File", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
trajectory_value = Tkinter.Entry(trajectory_frame, bd = 5)

wrapped = Tkinter.IntVar()
wrapped_frame = Tkinter.Frame(top)
wrapped_label = Tkinter.Label(wrapped_frame, text = "", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
wrapped_value = Tkinter.Checkbutton(wrapped_frame, text = "wrapped", justify = Tkinter.LEFT, variable = wrapped)

startframe_frame = Tkinter.Frame(top)
startframe_label = Tkinter.Label(startframe_frame, text = "Start Frame", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
startframe_value = Tkinter.Spinbox(startframe_frame, from_=0, to_=10**100)

endframe_frame = Tkinter.Frame(top)
endframe_label = Tkinter.Label(endframe_frame, text = "End Frame", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
endframe_value = Tkinter.Spinbox(endframe_frame, from_=0, to_=10**100)

frameaverage_frame = Tkinter.Frame(top)
frameaverage_label = Tkinter.Label(frameaverage_frame, text = "Number of Frames\nto Average", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
frameaverage_value = Tkinter.Spinbox(frameaverage_frame, from_=1, to_ = 10**100)

timestep_frame = Tkinter.Frame(top)
timestep_label = Tkinter.Label(timestep_frame, text = "Trajectory Time Step", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
timestep_value = Tkinter.Entry(timestep_frame, bd = 5)

gro_frame = Tkinter.Frame(top)
gro_label = Tkinter.Label(gro_frame, text = "Gro File (optional)", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
gro_value = Tkinter.Entry(gro_frame, bd = 5)

dimension_frame = Tkinter.Frame(top)
dimension_label = Tkinter.Label(dimension_frame, text = "Dimension", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
dimension_value = Tkinter.Spinbox(dimension_frame, from_=1, to_ = 100)
dimension_value.delete(0,"end")
dimension_value.insert(0,3)

### Computation Options
atomtypes_frame = Tkinter.Frame(top)
atomtypes_label = Tkinter.Label(atomtypes_frame, text = "Atom Types", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
atomtypes_value = Tkinter.Entry(atomtypes_frame, bd = 5)

atomtypes2_frame = Tkinter.Frame(top)
atomtypes2_label = Tkinter.Label(atomtypes2_frame, text = "Atom Types 2", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
atomtypes2_value = Tkinter.Entry(atomtypes2_frame, bd = 5)

scatteringlengths_frame = Tkinter.Frame(top)
scatteringlengths_label = Tkinter.Label(scatteringlengths_frame, text = "Scattering Lengths", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
scatteringlengths_value = Tkinter.Entry(scatteringlengths_frame, bd = 5)

scatteringlengths2_frame = Tkinter.Frame(top)
scatteringlengths2_label = Tkinter.Label(scatteringlengths2_frame, text = "Scattering Lengths 2", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
scatteringlengths2_value = Tkinter.Entry(scatteringlengths2_frame, bd = 5)

atomgroup_frame = Tkinter.Frame(top)
atomgroup_label = Tkinter.Label(atomgroup_frame, text = "Atom Group", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
atomgroup_value = Tkinter.Entry(atomgroup_frame, bd = 5)

atomgroup2_frame = Tkinter.Frame(top)
atomgroup2_label = Tkinter.Label(atomgroup2_frame, text = "Atom Group 2", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
atomgroup2_value = Tkinter.Entry(atomgroup2_frame, bd = 5)

numberbins_frame = Tkinter.Frame(top)
numberbins_label = Tkinter.Label(numberbins_frame, text = "Number of Bins", width = 20, justify = Tkinter.LEFT, anchor=Tkinter.W)
numberbins_value = Tkinter.Spinbox(numberbins_frame, from_=1, to_ = 10**10)
numberbins_value.delete(0,"end")
numberbins_value.insert(0,100)

cutofflength_frame = Tkinter.Frame(top)
cutofflength_label = Tkinter.Label(cutofflength_frame, text = "Max Cutoff Length", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
cutofflength_value = Tkinter.Entry(cutofflength_frame, bd = 5)

kstart_frame = Tkinter.Frame(top)
kstart_label = Tkinter.Label(kstart_frame, text = "k Start Index" , width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
kstart_value = Tkinter.Spinbox(kstart_frame, from_=0, to_ = 10**10)

timepoints_frame = Tkinter.Frame(top)
timepoints_label = Tkinter.Label(timepoints_frame, text = "Number of Output\nTime Points", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
timepoints_value = Tkinter.Spinbox(timepoints_frame, from_ = 1, to_ = 10**10)
timepoints_value.delete(0,"end")
timepoints_value.insert(0,100)

timescale = Tkinter.StringVar(top)
timescale.set("linear")
timescale_frame = Tkinter.Frame(top)
timescale_label = Tkinter.Label(timescale_frame, text = "Output Timescale", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
timescale_value = Tkinter.OptionMenu(timescale_frame, timescale, "linear", "log")

interval_frame = Tkinter.Frame(top)
interval_label = Tkinter.Label(interval_frame, text = "Frame Interval" , width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
interval_value = Tkinter.Entry(interval_frame, bd = 5)

overlap_frame = Tkinter.Frame(top)
overlap_label = Tkinter.Label(overlap_frame, text = "Overlap Length", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
overlap_value = Tkinter.Entry(overlap_frame, bd = 5)

is_averaged = Tkinter.IntVar()
qaverage_frame = Tkinter.Frame(top)
qaverage_label = Tkinter.Label(qaverage_frame, text = "", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
qaverage_value = Tkinter.Checkbutton(qaverage_frame, text = "q average", justify = Tkinter.LEFT, variable = is_averaged)

add_bar = Tkinter.IntVar()
qbar_frame = Tkinter.Frame(top)
qbar_label = Tkinter.Label(qbar_frame, text = "", width = 20, justify = Tkinter.LEFT, anchor = Tkinter.W)
qbar_value = Tkinter.Checkbutton(qbar_frame, text = "q bar", justify = Tkinter.LEFT, variable = add_bar)

bondorder_frame = Tkinter.Frame(top)
bondorder_label = Tkinter.Label(bondorder_frame, text = "Bond Order", width = 20, justify =  Tkinter.LEFT, anchor = Tkinter.W)
bondorder_value = Tkinter.Spinbox(bondorder_frame, from_ = 1, to_ = 10**10)
bondorder_value.delete(0,"end")
bondorder_value.insert(0,6)

### Base Menu
quit_frame = Tkinter.Frame(top)
quit_frame.pack(in_=top, side= Tkinter.BOTTOM)
quit_button = Tkinter.Button(quit_frame, text = "Quit", fg = "red",  command = top.quit)
quit_button.pack(in_ = quit_frame, side = Tkinter.BOTTOM)
submit_frame = Tkinter.Frame(top)
submit_frame.pack(in_=top, side= Tkinter.BOTTOM, pady=15)
input_button = Tkinter.Button(submit_frame, text = "Generate Input Script", fg = "blue",  command = generate_input_file)
input_button.pack(in_=submit_frame, side = Tkinter.LEFT, padx = 5)
submit_button = Tkinter.Button(submit_frame, text = "Compute Quantity", fg = "blue",  command = compute_quantity)
submit_button.pack(in_=submit_frame, side = Tkinter.LEFT, padx = 5)

### operation
top.mainloop()
top.destroy()