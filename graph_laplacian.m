
function L = graph_laplacian(W)

DCol = full(sum(W,2));
D = spdiags(DCol,0,speye(size(W,1)));
L = D - W;