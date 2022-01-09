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

# Environment setup
 
To install, first follow the readme instructions inside ompbox1, ompbox10, and ksvdbox13. Certain of their elements must be compiled locally. This is done by entering, for example, /ksvdbox13/private and running "matlab make". See the readme inside each folder for details.

Then, ensure that those paths are exposed to your MATLAB environment using 'addpath'. This incldues the GraphSC and DGRDL_Package folder as well.

Finally, we've found that sometimes, functions in DGRDL_Package still can't see functions that it needs in ksvdbox13, namely the sprow function. One workaround is, post-compilation, copy ksvdbox13/private/sprow.* to DGRDL_Package/. 

# Helper scripts

The helper scripts, DGRDL_wrapper and OMP_wrapper, are simple wrappers that I wrote in order to facilitate development and scripting for the Webster paper. I incorporated the graph-laplacian calculating script from GraphSC, so that given an input data matrix, those calculations are performed automatically and passed on to DGRDL. I also added a k-medioids step in dictionary initialization, such that any DGRDL run begins with an representative and stable set of k elements as the initial dictionary (rather than randomized columns, which was the default in the original implementation). Finally, these wrappers obey the data conventions established in the k-SVD, OMP, DGRDL and GSC implementations. For your input data matrix:

1. Features are rows
2. Signals are columns.

For image analysis, this means that your data matrix must be in the form pixels x images. For Webster, your data matrix must be in the form cells x genes.

Finally, each of these methods expects a data matrix as input, so make sure when you import your matrix into Matlab, that it has the right format and that column and row names are stored appropriately. For DGRDL_wrapper, I found it easiest to pass along the path to a *completely numerical* matrix without rownames and column names. This worked for me since I was using R to export flat files, read back the outputs, and put the row and column names back in through my R session. Your mileage may vary.

The output of both functions is a .mat object that contains all the inputs and outputs of the runs.

# Note on MATLAB licenses
MATLAB requires a license to run, so contact your institutional IT department if you need help obtaining a license for running these scripts. For those without access, there is a free trial (https://www.mathworks.com/campaigns/products/trials.html). We are also working on a Docker Image that contains a compiled version of DGRDL_wrapper and OMP_wrapper, for those who wish to run those scripts specifically. 
