function [D,X,L] = DGRDL_optL(params)
%=============================================
% Dual graph regularized dictionary learning algorithm (DGRDL)
% with optimized Laplacian.
%
% DGRDL_optL solves the optimization problem:
%       min  |Y-D*X|_F^2 + alpha*tr(D'*L*D) + beta*tr(X*Lc*X') + mu*|L|_F^2  
%      D,X,L
%       s.t.  |x_i|_0 <= T
%             L(i,j)=L(j,i)<=0 (i~=j)
%             L*ones(N,1) = zeros(N,1)
%             trace(L) = N
%
% DGRDL parameters:
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
% Sparse coding parameters:
%       SCiternum - number of ADMM iterations
%       rho - ADMM step size parameter
%    
% Graph learning parameters:
%       mu - regularization coefficient
%       Liter - frequency of Laplacian updates (default: 4)
%       thr - weight threshold (discard edges with weight<thr)
%
% Outputs:
%       D - learned dictionary
%       X - sparse coefficient matrix
%       L - learned graph Laplacian
%
% References:
% Y. Yankelevsky and M. Elad, "Dual Graph Regularized Dictionary Learning",  
% IEEE Transactions on Signal and Information Processing over Networks, 
% vol. 2, no. 4, pp. 611-624, Dec. 2016.
%
% Yael Yankelevsky
% Computer Science Department
% Technion - IIT
% yaelyan@cs.technion.ac.il
% June 2016
%=============================================

if ~isfield(params,'L')
    error('Initial graph Laplacian L missing!');
end

if isfield(params,'alpha')
    alpha = params.alpha;
else
    error('Regularizaion coefficient alpha missing!');
end

if isfield(params,'mu')
    mu = params.mu;
else
    error('Regularizaion coefficient mu missing!');
end

if isfield(params,'iternum')
    iternum = params.iternum;
else
    iternum = 20;
end

if isfield(params,'Liter')
    Liter = params.Liter;
else
    Liter = 4;
end

innerIter = round(iternum/Liter);
params.iternum = innerIter;
Lparams.alpha = alpha;
Lparams.mu = mu;
for iter=1:Liter
    fprintf('Global update iteration %d/%d\n',iter,Liter);
    [curD,X] = DGRDL(params);
    curD = normcols(curD);
    Lparams.D = curD;
    params.D = curD;
    curL = optimizeGraphLaplacian(Lparams);
    params.L = curL;
end
D = curD;
L = curL;
end
