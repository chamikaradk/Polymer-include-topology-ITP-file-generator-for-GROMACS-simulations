# Polymer-ITP-generator
Generate automated initial topology for gromacs simulations

This code automatically generate polymer itp file for gromacs simulations. Monomer itp file needs to be input, along with bonded parameters of 
linkers, then assign repeat units to generate the polymer itp. We rely on the atom indexing in the following way - Nn = N1+n x N+1H
Materials studio can be used to generate the polymeric structure from monomer pdb.
Use lig Par gen server by jorgensen lab @ Yale to generate initial itp file for the monomer

Input the following information in the code to generate your polymer ITP combining explicit calculations of dihedrals and charges
