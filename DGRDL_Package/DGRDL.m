function [D,X] = DGRDL(params)
%=============================================
% Dual graph regularized dictionary learning algorithm (DGRDL).
%
% DGRDL solves the optimization problem:
%       min  |Y-D*X|_F^2 + alpha*tr(D'*L*D) + beta*tr(X*Lc*X')  s.t.  |x_i|_0 <= T
%       D,X
%
% Main parameters:
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
% Outputs:
%       D - learned dictionary
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

if isfield(params,'Y')
    Y = params.Y; 
else
    error('Input data matrix Y missing!');
end
M = size(Y,2); % number of signals

if (isfield(params,'D'))
    D = params.D;
    K = size(D,2); % number of atoms
    if (isfield(params,'K'))
        if (params.K ~= K)
            error('Invalid initial dictionary D: size does not match the field K');
        end
    end
else
    if (isfield(params,'K'))
        K = params.K;
        % initialize the dictionary
        data_ids = find(colnorms_squared(Y)>1e-6);   % ensure no zero data elements are chosen
        perm = randperm(length(data_ids));
        D = Y(:,data_ids(perm(1:K)));
    else
        error('Missing input dictionary or dictionary size!');
    end
end
D = normcols(D); % make sure the atoms are normalized

if isfield(params,'T')
    T = params.T;
else
    error('Sparsity constraint T missing!');
end

if isfield(params,'L')
    L = params.L;
else
    error('Graph Laplacian L missing!');
end

if isfield(params,'Lc')
    Lc = params.Lc;
else
    error('Manifold Laplacian Lc missing!');
end

if isfield(params,'alpha')
    alpha = params.alpha;
else
    error('Regularizaion coefficient alpha missing!');
end

if isfield(params,'beta')
    beta = params.beta;
else
    error('Regularizaion coefficient beta missing!');
end

if isfield(params,'iternum')
    iternum = params.iternum;
else
    iternum = 20;
end

if isfield(params,'sdebug') % activate debug mode
    sdebug = params.sdebug;
else
    sdebug = 0;
end

%TODO: Add dimension check on Laplacian


% sparse coding parameters
scparams.D = D;
scparams.Y = Y;
scparams.T = T;
scparams.L = Lc;
scparams.beta = beta;
if isfield(params,'SCiternum')
    scparams.iternum = params.SCiternum;
end
if isfield(params,'rho')
    scparams.rho = params.rho;
end
muthresh = 0.99; 

for iter=1:iternum
    fprintf('Dictionary update iteration %d/%d\n',iter,iternum);
    % sparse coding
    X = GRSC_ADMM(scparams);
    % dictionary update
    replaced_atoms = zeros(1,K);  % mark each atom replaced by optimize_atom.
    unused_sigs = 1:M;  % tracks the signals that were used to replace "dead" atoms. 
                        % makes sure the same signal is not selected twice.
    atomperm = randperm(K);
    for j=1:K
        [D(:,atomperm(j)),x_j,data_indices,unused_sigs,replaced_atoms] = ...
            optimize_atom(Y,D,atomperm(j),X,unused_sigs,replaced_atoms,L,Lc,alpha,beta);
        X(atomperm(j),data_indices) = x_j; % plug back into X
    end
    
    [D,cleared_atoms] = cleardict(D,X,Y,muthresh,unused_sigs,replaced_atoms);
    totalReplaced = sum(replaced_atoms) + cleared_atoms;
    if (totalReplaced>0 && sdebug)
        fprintf('   replaced %d atoms\n', totalReplaced);
    end
    
    scparams.D = D;
end
end

     
function [atom,x_j,data_indices,unused_sigs,replaced_atoms] = ...
    optimize_atom(Y,D,j,X,unused_sigs,replaced_atoms,L,Lc,alpha,beta)
[x_j,data_indices] = sprow(X,j); 
smallX = X(:,data_indices);
Dj = D(:,j);

if (length(data_indices) < 1) %found unused atom
    maxsignals = 5000;
    perm = randperm(length(unused_sigs));
    perm = perm(1:min(maxsignals,end));
    E = sum((Y(:,unused_sigs(perm)) - D*X(:,unused_sigs(perm))).^2);
    [~,i] = max(E);
    atom = Y(:,unused_sigs(perm(i)));
    atom = atom./norm(atom);
    x_j = zeros(size(x_j));
    unused_sigs = unused_sigs([1:perm(i)-1,perm(i)+1:end]);
    replaced_atoms(j) = 1;
    return;
end

if ~isempty(data_indices)
    EjR = Y(:,data_indices)-D*smallX+Dj*x_j;
    % update atom j
    atom = (norm(x_j)^2*eye(size(L))+alpha*L)\(EjR*x_j');
    atom = atom/norm(atom);
    % update x_j
    x_j = (eye(length(data_indices))+beta*Lc(data_indices,data_indices))\(EjR'*atom);
    x_j = x_j';
else
    atom = Dj;
end
end


function Y = colnorms_squared(X)
% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
    blockids = i : min(i+blocksize-1,size(X,2));
    Y(blockids) = sum(X(:,blockids).^2);
end
end


function [D,cleared_atoms] = cleardict(D,X,Y,muthresh,unused_sigs,replaced_atoms)
use_thresh = 4;  % at least this number of samples must use the atom to be kept
K = size(D,2);

% compute error in blocks to conserve memory
err = zeros(1,size(Y,2));
blocks = [1:3000:size(Y,2) size(Y,2)+1];
for i = 1:length(blocks)-1
    err(blocks(i):blocks(i+1)-1) = sum((Y(:,blocks(i):blocks(i+1)-1)-D*X(:,blocks(i):blocks(i+1)-1)).^2);
end

cleared_atoms = 0;
usecount = sum(abs(X)>1e-7, 2);
for j = 1:K
    % compute G(:,j)
    Gj = D'*D(:,j);
    Gj(j) = 0;
    % replace atom
    if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
        [~,i] = max(err(unused_sigs));
        D(:,j) = Y(:,unused_sigs(i)) / norm(Y(:,unused_sigs(i)));
        unused_sigs = unused_sigs([1:i-1,i+1:end]);
        cleared_atoms = cleared_atoms+1;
    end
end
end
