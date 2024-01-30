# Polymer-ITP-generator (ITP_Gen)
PolyITPgen Generates automated initial polymer topology for gromacs simulations. Both periodic and non periodic polymer molecules can be generated.
Limitations: Only work with monomers with two interaction sites, makig linear polymer chains. 
This code automatically generate polymer itp file for gromacs simulations from Monomer itp file. You may use LigParGen server by jorgensen lab @ Yale to generate initial itp file for the monomer.
Monomer needs to be input, along with bonded interaction parameters pertaining to:
  1.Inter molecular and; 
  2.Intra molecular properties (linkers). 
Then, simply assign number of repeat units to generate the polymer itp. The code relies on the atom indexing to start and end with the two terminal H atoms of the monomer. 
PolGroGen can be used to generate the polymeric structure input file from monomer pdb - Please check the other repositary PoGroGen.

