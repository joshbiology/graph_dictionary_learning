# Graph Dictionary Learning
This repository is a convenience repo for anyone interested in applying two graph regularized dictionary learning strategies:
1. Dual Graph Regularized Dictionary Learning (DGRDL), a graph regularized dictionary learning method from Yael Yankelevsky and Michael Elad. 10.1109/TSIPN.2016.2605763
2. Graph Regularized Sparse Coding (GSC), a graph regularized dictionary learning method from Miao Zheng, Deng Cai et al. 10.1109/TIP.2010.2090535

DGRDL is the core factorization method underlying Webster, an approach for sparsely factorizing fitness matrices. More on that here: https://github.com/joshbiology/gene_fn. 

DGRDL is built on top of the classic dictionary learning approach, k-SVD (k-Singular Value Decomposition) and OMP (Orthogonal Matching Pursuit). The specific implementations are from this paper:

3. Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit, (http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf, direct download of the paper).

The source code for those papers was obtained from these websites:

1. DGRDL: https://elad.cs.technion.ac.il/software/
2. GSC: http://www.cad.zju.edu.cn/home/dengcai/Data/SparseCoding.html
3. K-SVD/OMP: http://www.cs.technion.ac.il/~ronrubin/software.html

# Overview
 
To install, first follow the readme instructions inside ompbox1, ompbox10, and ksvdbox13. Certain of their elements must be compiled locally. Then, ensure that those paths are exposed to your MATLAB environment using 'addpath'. 

The helper scripts, DGRDL_wrapper and OMP_wrapper, are simple wrappers that I wrote in order to facilitate development and scripting for the Webster paper. I incorporated the graph-laplacian calculating script from GSC, so that given an input data matrix, those calculations are performed automatically and passed on to DGRDL. I also added a k-medioids step in dictionary initialization, such that any DGRDL run begins with an representative and stable set of k elements as the initial dictionary (rather than randomized columns, which is the default). Finally, these wrappers obey the convention established in the k-SVD, OMP, DGRDL and GSC implementations which is that for your input data matrix:

1. Features are rows
2. Signals are columns.

For image analysis, this means that your data matrix must be in the form pixels x images. For Webster, your data matrix must be in the form cells x genes.

Finally, each of these methods expects a data matrix as input, so make sure when you import your matrix into Matlab, that it has the right format and that column and row names are stored appropriately. For DGRDL_wrapper, I found it easiest to pass along the path to a *completely numerical* matrix without rownames and column names. This worked for me since I was using R to export flat files, read back the outputs, and put the row and column names back in through my R session. Your mileage may vary.