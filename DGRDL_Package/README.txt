README for DGRDL Matlab Package Version 1.1.

This is the Matlab Package for the Dual Graph Regularized Dictionary Learning (DGRDL)
algorithm, presented in:
Y. Yankelevsky and M. Elad, "Dual Graph Regularized Dictionary Learning",  
IEEE Transactions on Signal and Information Processing over Networks, 
vol. 2, no. 4, pp. 611-624, Dec. 2016.

INSTALLATION:
- The DGRDL code requires to have the OMP-BOX and KSVD-BOX installed, which 
can be freely downloaded from
http://www.cs.technion.ac.il/~ronrubin/software.html
- Once these packages have their corresponding mex files compiled 
(check their instructions for further details), make sure they are added in
the Matlab path.

USAGE:
The main function is DGRDL.m, which performs dictionary learning on the 
indicated training data and outputs a graph dictionary (refer to the 
referenced paper for more details). 
The function GRSC_ADMM.m computes the manifold regularized sparse codes 
(see algorithm 2 in the referenced paper).
To jointly optimize the graph with the dictionary, use the function DGRDL_optL.m.

All comments are welcome at yaelyan@cs.technion.ac.il

Yael Yankelevsky
Computer Science Department - Technion
June, 2016.