function L = optimizeGraphLaplacian(params)
%=============================================
% Optimize the graph Laplacian matrix based on dictionary atoms regularity
%
% The function solves the optimization problem:
%   min  alpha*trace(Y'*L*Y)+mu*|L|_F^2   s.t.  L(i,j)=L(j,i)<=0 (i~=j)
%    L                                          L*ones(N,1) = zeros(N,1)
%                                               trace(L) = N
%
% Main parameters:
%       D - input dictionary/data matrix
%       alpha,mu - regularization coefficients (defaults: 1,0.8)
%       thr - weight threshold - discard edges with weight<thr (default: 0)
%       useCVX - optimization mode (default: 0)
%    
% Output:
%       L - NxN Laplacian matrix
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

if isfield(params,'D')
    D = params.D;
else
    error('Input D missing!');
end

if (isfield(params,'alpha'))
    alpha = params.alpha;
else
    alpha = 1;
    params.alpha = alpha;
end
if (isfield(params,'mu'))
    mu = params.mu;
else
    mu = 0.8;
    params.mu = mu;
end

if (isfield(params,'thr'))
    thr = params.thr;
else
    thr = 0;
    params.thr = thr;
end

if (isfield(params,'useCVX'))
    useCVX = params.useCVX;
else
    useCVX = 0;
end


N = size(D,1);

if (useCVX)
    if ~exist('cvx_begin')
        error('cvx Package missing!');
    end
    cvx_begin
        variable L(N,N) symmetric;
        minimize(alpha*trace(D'*L*D)+mu*square_pos(norm(L,'fro')));
        subject to
            L*ones(N,1) == zeros(N,1);
            L-diag(diag(L))<=0;
            trace(L) == N;
    cvx_end
else    
    [I,J] = find(tril(ones(N)));
    loctri = find(I~=J & J<=N & I<=N);
    It = sub2ind([N,N], J(loctri), I(loctri));
    I = sub2ind([N,N],I,J);
    M = sparse([I;It],[(1:length(I))';loctri],1,N^2,length(I));
    
    A = sparse(N,N*(N+1)/2);
    Sk=0;
    for i=1:N
        Si = N-i+1;
        A(i,Sk+1:Sk+Si) = 1;
        A(i+1:N,Sk+2:Sk+Si) = eye(N-i);
        Sk = Sk+Si;
    end
    diagB = ones(N*(N+1)/2,1);
    Sk=0;
    for k=N:-1:1
        diagB(Sk+1)=0;
        Sk=Sk+k;
    end
    B = spdiags(diagB,0,N*(N+1)/2,N*(N+1)/2);
    
    C = ones(1,N*(N+1)/2)*(speye(N*(N+1)/2)-B);
    ddvec = reshape(D*D',N^2,1);
    x = quadprog(2*mu*(M'*M),M'*ddvec*alpha,B,zeros(N*(N+1)/2,1),[A;C],[zeros(N,1);N]);
    L = zeros(N);
    L(tril(true(N))) = x;
    L = L+L'-diag(diag(L));
end
% discard weak weights
if (thr>0) 
    L(abs(L)<thr*max(abs(L(:)))) = 0;
    L(eye(N)==1)=(sum(L,2)-diag(L));
end
end