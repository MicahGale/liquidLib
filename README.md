# LiquidLib

## Purpose
A comprehensive toolbox for analyzing classical and ab initio molecular dynamics simulations of liquids and liquid-like matter with applications to neutron scattering experiments

## Introduction
Neutron scattering is a powerful experimental technique for characterizing the structure and dynamics of materials on the atomic or molecular scale. However, the interpretation of experimental data from neutron scattering is oftentimes not trivial, partly because scattering methods probe ensemble-averaged information in the reciprocal space. Therefore, computer simulations, such as classical and \textit{ab initio} molecular dynamics, are frequently used to unravel the time-dependent atomistic configurations that can reproduce the scattering patterns and thus assist in the understanding of the microscopic origin of certain properties of materials. LiquidLib is a post-processing package for analyzing the trajectory of atomistic simulations of liquids and liquid-like matter with application to neutron scattering experiments. From an atomistic simulation, LiquidLib provides the computation of various statistical quantities including the pair distribution function, the weighted and unweighted structure factor, the mean squared displacement, the non-Gaussian parameter, the four-point correlation function, the velocity auto correlation function, the self and collective van Hove correlation function, the self and collective intermediate scattering function, and the bond orientational order parameter. LiquidLib analyzes atomistic trajectories generated from packages such as LAMMPS, GROMACS, and VASP. It also offers an extendable platform to conveniently integrate new quantities into the library and integrate simulation trajectories of other file formats for analysis. Weighting the quantities by element-specific neutron-scattering lengths provides results directly comparable to neutron scattering measurements. Lastly, LiquidLib is independent of dimensionality, which allows analysis of trajectories in two, three, and higher dimensions. The code is beginning to find worldwide use.

## Capabilities
LiquidLib provides computation the following quantities:

* Pair Distribution Function
* Weighted and Unweighted Structure Factor
* Mean Squared Displacement
* Non-Gaussian Parameter
* Four-point Correlation Function
* Velocity Auto Correlation Function
* Self van Hove Correlation Function
* Collective van Hove Correlation Function
* Self Intermediate Scattering Function
* Collective Intermediate Scattering Function
* Bond Orientational Order Parameter
* And many more to come

And provides the further quantities in developed post publication:

* Center of Mass Velocity Auto Correlation Function
* Current Correlation Function
* Electric Current Correlation Function
* Trajectory Type Converter

Since the publication of the original LiquidLib, we have also added options for negative k vectors in k space quantities, progressive trajectory reading for memory limited computers or extremely large trajectories, a gui for easy input script construction, and fixed several bugs. These additions will be published in an extension paper soon and can currently be found in the develop branch on github.

## Reference
if used, please cite the following reference:

* N. P. Walter, A. Jaiswal, Z. Cai, and Y. Zhang, “LiquidLib: A comprehensive toolbox for analyzing classical and ab initio molecular dynamics simulations of liquids and liquid-like matter with applications to neutron scattering experiments,” Comput. Phys. Commun., vol. 228, pp. 209–218, 2018.

## Contact
This project is developed and maintained by Z Lab. Developers to this project were made by:

* Nathan Walter (@walternathan6754)
* Zhikun Cai (@caizkun)
* Zhixia Li (@zhixia721)
* Yanqin Zhai (@yanqinz2)

Previous developers
* Abhishek Jaiswal (@jaisabhi)

For obtaining user support please contact us at Z Lab.
