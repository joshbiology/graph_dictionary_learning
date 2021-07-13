function DGRDL_wrapper(filename_Y, varargin)
%=============================================
% Wrapper for Dual graph regularized dictionary learning algorithm (DGRDL).
% Convenience function for Pan et al., 2021
%
%       D - initial dictionary
%       K - dictionary size (if specified without the paramterer D, 
%                            the initial dictionary would be K randomly
%                            selected training signals)
%       Y - training data
%       T - sparsity constraint (max. number of coefficients per signal)
%       L - graph Laplacian
%       Lc - manifold graph Laplacian
%       alpha,beta - regularization coefficients
%       iternum - number of iterations (default: 20)
%    
%    
% Outputs:
%       D - learned dictionary
%       X - sparse coefficient matrix
%
% References:
% Y. Yankelevsky and M. Elad, "Dual Graph Regularized Dictionary Learning",  
% IEEE Transactions on Signal and Information Processing over Networks, 
% vol. 2, no. 4, pp. 611-624, Dec. 2016.


%=============== Argument validation ====================
default_filename_out = "out"; 
default_K = 300;
default_T = 4;
default_alpha = 0.2; 
default_beta = 0.6;
default_iternum = 20; 
default_seed = 1;
default_num_neighbor_gene = 5;
default_num_neighbor_cl = 5;
default_knn_weight_mode = 'Cosine';

p = inputParser;
validPosNum= @(x) isnumeric(x) && (x >= 0);
addRequired(p,'filename_Y',@ischar);
addParameter(p,'filename_out',default_filename_out,@ischar);
addParameter(p,'K',default_K,@isnumeric);
addParameter(p,'T',default_T,@isnumeric);
addParameter(p,'alpha',default_alpha,validPosNum);
addParameter(p,'beta',default_beta,validPosNum);
addParameter(p,'iternum',default_iternum,@isnumeric);
addParameter(p,'seed',default_seed,@isnumeric);
addParameter(p,'filename_init_D',"",@ischar);
addParameter(p,'num_neighbor_cl',default_num_neighbor_cl,@isnumeric);
addParameter(p,'num_neighbor_gene',default_num_neighbor_gene,@isnumeric);
addParameter(p,'knn_weight_mode',default_knn_weight_mode,@ischar);


parse(p,filename_Y,varargin{:});


%=============== Parameter set up====================

%Set up fixed seed for k-medioids. This allows reproducible output between
%runs of the same choice of K.
rng('default');

%Read in data; expect genes as columns.
fea = readmatrix(p.Results.filename_Y);

params = struct();
params.Y = fea;
params.T = p.Results.T;
params.K = p.Results.K;
params.alpha = p.Results.alpha;
params.beta = p.Results.beta;
params.iternum = p.Results.iternum;


%Set initial dictionary
if (p.Results.filename_init_D ~= "")
    fprintf('Initialized dict using %s\n',p.Results.filename_init_D);
    medoid_table = readtable(p.Results.filename_init_D);
    
    if (height(medoid_table)<p.Results.K)
        error('Incorrect dimension in medoid table.');
    else
        medoids = medoid_table{:,2}(1:p.Results.K);
    end
else
    fprintf('Initialized dict using k-medioids.\n');
    [~,~,~,~,medoids] = kmedoids(params.Y', p.Results.K, 'Distance', 'cosine');
end

params.D = fea(:,medoids);

% signal-signal graph (genes)
options = [];
options.WeightMode = p.Results.knn_weight_mode;
options.k = p.Results.num_neighbor_gene;

params.Lc = graph_laplacian(constructW(fea', options));

% internal graph (cell_lines)
options2 = [];
options2.WeightMode = p.Results.knn_weight_mode;
options2.k = p.Results.num_neighbor_cl;

params.L = graph_laplacian(constructW(fea, options2));

%=============== Run and save ====================

%Set up random seed for DGRDL. This allows seed to influence output of
%factorization but not centroid choice above.
rng(p.Results.seed);

[D,X] = DGRDL(params);

out = struct('D', D, 'X', X, 'K', p.Results.K, 'medoids', medoids,...
    'T', p.Results.T, 'alpha', p.Results.alpha, 'beta', p.Results.beta,...
    'iternum', p.Results.iternum, 'seed', p.Results.seed,'num_neighbor_cl', p.Results.num_neighbor_cl, ...
    'num_neighbor_gene', p.Results.num_neighbor_gene, 'filename_init_D', p.Results.filename_init_D,  ...
    'knn_weight_mode', p.Results.knn_weight_mode,...
    'filename_Y', p.Results.filename_Y, 'filename_out', p.Results.filename_Y);
save(p.Results.filename_out,'-struct', 'out');
end