# Displacement-controlled Newton-Raphson algorithm

## Description
This code calculates the nonlinear response (force-displacement curve) assuming material nonlinearity and geometric linearity using the displacement-controlled Newton-Raphson algorithm.
The input files are (nodes.txt, restraints.txt, elements_information.txt) defining the geometry of the system, the restraints of the nodes and the properties of the elements. The outputs
of main.m are .txt files for the displacement U, the corresponding force F and the residual force R for all iterations. The code plot_force_displacement.m creates a graph showing the
nonlinear response of the structure (i.e. U and F generated in main.m), the code plot_residual.m creates a graph showing the residual force for all iterations over all loadsteps
(i.e. R generated in main.m).

## Getting Started
The code was developped on Matlab version R2021b.

### Executing the program
* Make sure the input files (nodes.txt, restraints.txt, elements_information.txt) are in the same folder as main.m (otherwise modify the paths in the beginning of the code).
* Run "main.m" (creates output files U.txt, F.txt and R.txt)
* To plot the force-displacement curves: run plot_force_displacement.m (indicated the node to be plotted)
* To plot the residual force: run plot_residual.m

## Authors
[Nicole Widmer](nicole.widmer@epfl.ch)

##License
The project is licensed under the CC BY 4.0 license.

## Acknowledgments
The function get_nodal_forces is based on a code provided in the course "CIVIL-449 Nonlinear analysis of structures" at EPFL taught by the professors Katrin Beyer and Dimitrios Lignos.

 
