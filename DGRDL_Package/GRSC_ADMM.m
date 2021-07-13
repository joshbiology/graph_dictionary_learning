function X = GRSC_ADMM(params)
%=============================================
% Graph regularized sparse coding algorithm (GRSC_ADMM).
%
% GRSC_ADMM solves the optimization problem:
%       min  |Y-D*X|_F^2 + beta*tr(X*L*X')  s.t.  |x_i|_0 <= T
%        X
%
% Main parameters:
%       D - the dictionary (its columns MUST be normalized)
%       Y - the signals to represent
%       T - sparsity constraint (max. number of coefficients per signal)
%       L - manifold graph Laplacian
%       beta - regularization coefficient
%       X - initial sparse code (default: run non-regularized OMP)
%       iternum - number of ADMM iterations (default: 5)
%       rho - ADMM step size parameter (default: 1)
%       runDebiasing - update values on the determined support using least-squares (default: 1)
%   
% Output:
%       X - sparse coefficient matrix
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

if ~exist('omp','file')
    error('OMP Package missing!');
end

if isfield(params,'D')
    D = params.D;
    G = D'*D;
    if (sum(abs(diag(G)-1)>1e-12)>0)
        warning('Dictionary atoms are not normalized!');
    end
else
    error('Input dictionary D missing!');
end

if isfield(params,'Y')
    Y = params.Y;
else
    error('Input data matrix Y missing!');
end

if isfield(params,'T')
    T = params.T;
else
    error('Sparsity constraint T missing!');
end

if isfield(params,'L')
    L = params.L;
else
    error('Manifold Laplacian L missing!');
end

if isfield(params,'beta')
    beta = params.beta;
else
    error('Regularizaion coefficient beta missing!');
end

if isfield(params,'iternum')
    iternum = params.iternum;
else
    iternum = 5;
end

if isfield(params,'rho')
    rho = params.rho;
else
    rho = 1;
end

if (isfield(params,'runDebiasing'))
    runDebiasing = params.runDebiasing;
else
    runDebiasing = 1;
end
M = size(Y,2); % number of signals
K = size(D,2); % number of atoms

if isfield(params,'X')
    X = params.X;
else
    ompparams = {'checkdict','off'};
    X = omp(D'*Y,G,T,ompparams{:});
    if (beta==0)
        X = sparse(X);
        return
    end
end

Z = X;
U = zeros(K,M);

for i = 1:iternum
    X = sylvester(full(G+rho*eye(size(G))),full(beta*L),full(D'*Y+rho*(Z-U)));
    Z = SpProj(X+U,T);
    U = U+X-Z;
end
X = sparse(Z);

% update values on the determined support using LS
if (runDebiasing)
	for j=1:M
		suppInd = find(X(:,j)~=0);
		X(suppInd,j) = pinv(D(:,suppInd))*Y(:,j);
	end
end
end


function Z = SpProj(XU,T)
Z = zeros(size(XU));
for j=1:size(Z,2)
    [~,ind] = sort(abs(XU(:,j)),'descend');
    Z(ind(1:T),j) = XU(ind(1:T),j);
end
end
